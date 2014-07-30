#!/usr/bin/python
#
# This file is part of Mitsuba, a physically based rendering system.
#
# Copyright (c) 2007-2014 by Wenzel Jakob and others.
#
# Mitsuba is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License Version 3
# as published by the Free Software Foundation.
#
# Mitsuba is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from pprint import pprint
from boto.ec2.connection import EC2Connection
import os, sys, time, re, datetime, subprocess, base64, boto

# Configure this with the authentication data of your account
AWS_ACCESS_KEY_ID     = 'XXXXXXXXXXXXXXXXXXXX'
AWS_SECRET_ACCESS_KEY = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
KEYPAIR_NAME          = 'XXXXXXXXX'

# Where do you want to start nodes? (e.g. us-east-1, eu-west-1, etc..)
AWS_REGION            = 'us-east-1'

# Ensure that the SSH private key file of the following name is
# located in the same directory as this script
SSH_PRIVATE_KEY_FILE  = KEYPAIR_NAME + ".pem"

# Repository to be used for fetching Mitsuba Ubuntu packages
PKG_REPOSITORY        = 'deb https://www.mitsuba-renderer.org binary/'

# Optional: when syncing additional files to the cluster nodes,
# you can specify the relevant S3 bucket here
S3_PATH               = 's3://XXXXXX'

conn = None

# AMI ID: Stateless Ubuntu Maverick 64bit
ami_ids = {
	'ap-northeast-1' : 'ami-420fa443',
	'ap-southeast-1' : 'ami-12423c40',	
	'eu-west-1'	: 'ami-1b9ca96f',
	'us-east-1'	: 'ami-08f40561',
	'us-west-1'	: 'ami-a17e2ee4'
}

def remoteCommand(host, cmd, quiet = True):
	#print('Executing command "%s" on node %s' % (cmd, host))
	return subprocess.Popen(['ssh', '-i', SSH_PRIVATE_KEY_FILE, 
		'-o', 'StrictHostKeyChecking=no',
		'-o', 'LogLevel=ERROR',
		'-o', 'UserKnownHostsFile=/dev/null',
		'ubuntu@%s' % host,
		'bash -c \'%s\'' % cmd], stdout = subprocess.PIPE if quiet else sys.stdout)

def remoteAdminCommand(host, cmd, quiet = True):
	#print('Executing admin command "%s" on node %s' % (cmd, host))
	return subprocess.Popen(['ssh', '-i', SSH_PRIVATE_KEY_FILE, 
		'-o', 'StrictHostKeyChecking=no',
		'-o', 'LogLevel=ERROR',
		'-o', 'UserKnownHostsFile=/dev/null',
		'ubuntu@%s' % host,
		'sudo bash -c \'%s\'' % cmd], stdout = subprocess.PIPE if quiet else sys.stdout)

def parse_timestamp(s):
  """Returns (datetime, tz offset in minutes) or (None, None)."""
  m = re.match(""" ^
    (?P<year>-?[0-9]{4}) - (?P<month>[0-9]{2}) - (?P<day>[0-9]{2})
    T (?P<hour>[0-9]{2}) : (?P<minute>[0-9]{2}) : (?P<second>[0-9]{2})
    (?P<microsecond>\.[0-9]{1,6})?
    (?P<tz>
      Z | (?P<tz_hr>[-+][0-9]{2}) : (?P<tz_min>[0-9]{2})
    )?
    $ """, s, re.X)
  if m is not None:
    values = m.groupdict()
    if values["tz"] in ("Z", None):
      tz = 0
    else:
      tz = int(values["tz_hr"]) * 60 + int(values["tz_min"])
    if values["microsecond"] is None:
      values["microsecond"] = 0
    else:
      values["microsecond"] = values["microsecond"][1:]
      values["microsecond"] += "0" * (6 - len(values["microsecond"]))
    values = dict((k, int(v)) for k, v in values.iteritems()
                  if not k.startswith("tz"))
    try:
      return datetime.datetime(**values), tz
    except ValueError:
      pass
  return None, None
		
def addNodes(instanceType, nodeCount, groupName):
	print('Booting %i nodes of type %s (group name = "%s") ..' % (nodeCount, instanceType, groupName))
	ami_id = ami_ids[AWS_REGION]
	image = conn.get_image(ami_id)
	reservation = image.run(min_count=nodeCount, max_count=nodeCount,
			instance_type = instanceType, key_name = KEYPAIR_NAME,
			user_data=groupName)
	while True:
		time.sleep(2)
		runningCount = 0
		for i in reservation.instances:
			i.update()
			if i.state == u'running':
				runningCount += 1
		print("%i/%i nodes are ready." % (runningCount, nodeCount))
		if runningCount == nodeCount:
			print("All nodes are running. Note that the actual OS might still take a few minutes")
			print("until it is fully initialized. Please wait sufficiently long or the following")
			print("steps ('install', 'start', ..) will fail.")
			break;

def addSpotNodes(instanceType, nodeCount, maxPrice, groupName):
	print('Requesting %i spot nodes of type %s (group name = "%s", max. price=%f)..' % (nodeCount, instanceType, groupName, maxPrice))
	ami_id = ami_ids[AWS_REGION]
	conn.request_spot_instances(str(maxPrice), ami_id, count=nodeCount,
			instance_type = instanceType, key_name = KEYPAIR_NAME,
			user_data=groupName, type='one-time')
	print('Done.')

def getGroup(instance):
	return base64.b64decode(conn.get_instance_attribute(instance.id, 'userData')['userData'])

def status():
	print("Querying spot instance requests...")
	spotInstances = conn.get_all_spot_instance_requests()
	for i in spotInstances:
		print("  %s: status=%s, price must be <= %.3f$" % (i.id, i.state, i.price))
	print("")
	print("Querying instances ...")
	reservations = conn.get_all_instances()
	instances = [i for r in reservations for i in r.instances]
	nodesByGroup = {}
	print("")
	for i in instances:
		if i.state != 'terminated':
			group = getGroup(i)
			if not group in nodesByGroup:
				nodesByGroup[group] = []
			spot_str = ""
			if i.spot_instance_request_id != None:
				spot_str = ", spot request: %s" % i.spot_instance_request_id
			dt, tz = parse_timestamp(i.launch_time)
			curTime = datetime.datetime(*time.gmtime()[:6])
			delta = curTime-dt
			hours, remainder = divmod(delta.seconds, 3600)
			minutes, seconds = divmod(remainder, 60)
			nodesByGroup[group] += ["  %s is %s (type: %s, running for: %id %ih %im, internal IP: %s%s)" % (i.public_dns_name, i.state, i.instance_type,
				delta.days, hours, minutes, i.private_ip_address, spot_str)]
	for g in nodesByGroup:
		print('Nodes in group "%s"' % g)
		print('===================')
		for n in nodesByGroup[g]:
			print(n)

def terminate(name):
	reservations = conn.get_all_instances()
	instances = [i for r in reservations for i in r.instances]
	for i in instances:
		if i.state == u'running' and i.public_dns_name == name:
			print("Stopping node %s .." % i.public_dns_name)
			i.terminate()
			return
	print('Node could not be found or is not running')


def terminateAll(groupName):
	print('Terminating all nodes in group \"%s\" ..' % groupName)
	reservations = conn.get_all_instances()
	instances = [i for r in reservations for i in r.instances]
	for i in instances:
		if i.state == u'running' and getGroup(i) == groupName:
			print("Stopping node %s .." % i.public_dns_name)
			i.terminate()

def cancelSpot(id):
	print('Terminating spot request \"%s\" ..' % id)
	result = conn.cancel_spot_instance_requests([id])
	for n in result:
		if n.spot_instance_request_id == id:
			print('Success')
			return
	print('Failed!')

def cancelAllSpot():
	print('Terminating all spot requests ...')
	spotInstances = conn.get_all_spot_instance_requests()
	for i in spotInstances:
		if i.state == 'active' or i.state == 'open':
			cancelSpot(i.id)

def install(groupName):
	print('Installing Mitsuba on all nodes of group "%s" ..' % groupName)
	reservations = conn.get_all_instances()
	instances = [i for r in reservations for i in r.instances]
	processes = []
	for i in instances:
		if i.state == u'running' and getGroup(i) == groupName:
			print("Sending command to node %s" % i.public_dns_name)
			processes += [remoteAdminCommand(i.public_dns_name, 'echo \"%s\" > /etc/apt/sources.list.d/mitsuba.list; dpkg --purge mitsuba; apt-get clean; export DEBIAN_FRONTEND=noninteractive; apt-get update; apt-get -q -y --allow-unauthenticated install mitsuba s3cmd; chown ubuntu /mnt' % PKG_REPOSITORY)]
	while True:
		doneCount = 0
		time.sleep(2)
		for p in processes:
			if p.poll() != None:
				doneCount += 1
		print("%i/%i nodes are ready." % (doneCount, len(processes)))
		if doneCount == len(processes):
			print("All nodes are ready.")
			break;

def systemLoad(groupName):
	print('Querying system load on all nodes of group "%s" ..' % groupName)
	reservations = conn.get_all_instances()
	instances = [i for r in reservations for i in r.instances]
	processes = []
	for i in instances:
		if i.state == u'running' and getGroup(i) == groupName:
			processes += [remoteCommand(i.public_dns_name, 'echo `cat /proc/loadavg | cut --delimiter=" " --fields=1-3`  --  `ec2metadata --public-hostname`', False)]
	while True:
		doneCount = 0
		time.sleep(1)
		for p in processes:
			if p.poll() != None:
				doneCount += 1
		if doneCount == len(processes):
			print("Done.")
			break;

def runCommand(cmd, groupName):
	print('Executing command "%s" on all nodes of group "%s" ..' % (cmd, groupName))
	reservations = conn.get_all_instances()
	instances = [i for r in reservations for i in r.instances]
	processes = []
	for i in instances:
		if i.state == u'running' and getGroup(i) == groupName:
			print("Sending command to node %s" % i.public_dns_name)
			processes += [remoteCommand(i.public_dns_name, cmd, False)]
	while True:
		doneCount = 0
		time.sleep(1)
		for p in processes:
			if p.poll() != None:
				doneCount += 1
		if doneCount == len(processes):
			print("Done.")
			break;

def syncData(prefix, groupName):
	print('Fetching S3 prefix "%s/%s" on all nodes of the group "%s" ..' % (S3_PATH, prefix, groupName))
	reservations = conn.get_all_instances()
	instances = [i for r in reservations for i in r.instances]
	processes = []
	for i in instances:
		if i.state == u'running' and getGroup(i) == groupName:
			print("Sending command to node %s" % i.public_dns_name)
			os.system('scp -i %s -o UserKnownHostsFile=/dev/null -o LogLevel=ERROR -o StrictHostKeyChecking=no -q .s3cfg ubuntu@%s:' % (SSH_PRIVATE_KEY_FILE, i.public_dns_name))
			processes += [remoteCommand(i.public_dns_name, 'cd /mnt; s3cmd sync %s/%s .' % (S3_PATH, prefix))]
	while True:
		doneCount = 0
		time.sleep(2)
		for p in processes:
			if p.poll() != None:
				doneCount += 1
		print("%i/%i nodes are ready." % (doneCount, len(processes)))
		if doneCount == len(processes):
			print("All nodes are ready.")
			break;

def start(groupName):
	print('Creating a Mitsuba cluster using the nodes of group "%s"' % groupName)
	reservations = conn.get_all_instances()
	instances = [i for r in reservations for i in r.instances]
	activeNodes = []
	activeNodeIPs = []
	processes = []
	for i in instances:
		if i.state == u'running' and getGroup(i) == groupName:
			activeNodes += [i.public_dns_name]
			activeNodeIPs += [i.private_ip_address]
	if len(activeNodes) == 0:
		print("There are no running nodes!")
		return
	headNode = activeNodes[len(activeNodes)-1]
	list.remove(activeNodes, headNode)
	list.remove(activeNodeIPs, activeNodeIPs[len(activeNodeIPs)-1])
	for i in range(0, len(activeNodes)):
			print("Sending command to node %s" % activeNodes[i])
			processes += [remoteCommand(activeNodes[i], 'killall --quiet mtssrv; cd /mnt; nohup mtssrv -vq >/dev/null >&1 &')]
	connArgument = str.join(';', activeNodeIPs)
	if len(connArgument) > 0:
		connArgument = '-c \"' + connArgument + "\""
	while True:
		doneCount = 0
		time.sleep(2)
		for p in processes:
			if p.poll() != None:
				doneCount += 1
		print("%i/%i nodes are ready." % (doneCount, len(processes)))
		if doneCount == len(processes):
			print("All nodes are ready.")
			break;
	print("Creating head node ..")
	p = remoteCommand(headNode, 'killall --quiet mtssrv; cd /mnt; let pcount=`grep processor /proc/cpuinfo | wc -l`; let pcount=`perl -e "use POSIX qw/ceil/; print $pcount-2 < 1 ? 1 : $pcount-2"`; nohup mtssrv -p$pcount -vq %s >/dev/null >&1 &' % connArgument)
	p.wait()
	print('Done -- you can specify the head node "%s" in the Mitsuba network rendering dialog' % headNode)

def spotPrices(instanceType):
	start = datetime.datetime.now()
	end = datetime.datetime.now()
	start = start.replace(hour=0, minute=0, second=0, microsecond=0)
	end = end.replace(second=0, microsecond=0)

	startTime = start.isoformat() + str.format('Z{0:+06.2f}', float(time.timezone) / 3600)
	endTime = end.isoformat() + str.format('Z{0:+06.2f}', float(time.timezone) / 3600)
	history = conn.get_spot_price_history(start_time=startTime, end_time=endTime, 
			instance_type=instanceType, product_description = 'Linux/UNIX')

	for h in history:
		timestamp, tz = parse_timestamp(h.timestamp)
		print ("%s => %.5f dollars" % (timestamp, h.price))

def login(name):
	os.system('ssh -o UserKnownHostsFile=/dev/null -o LogLevel=ERROR -o StrictHostKeyChecking=no -i %s ubuntu@%s' % (SSH_PRIVATE_KEY_FILE, name))

if len(sys.argv) == 1:
	print('')
	print('                    Mitsuba Amazon EC2 instance launcher')
	print('')
	print(' Copyright (C) 2007-2014 by Wenzel Jakob and others. This is free software;')
	print(' see the source for copying conditions. There is NO warranty; not even for')
	print(' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.')
	print('')
	print(' Syntax: ./cluster.py <command> [arguments]    where command is one of')
	print('')
	print('   status -- Prints a list of all currently running cluster nodes')
	print('')
	print('   addNodes [instance-type] [number] <group> -- adds further machines to the')
	print('                pool of rendering nodes. The following instances types are')
	print('                curently available (N. Virginia pricing as of Jan 8 2011)')
	print('')
	print('                t1.micro    ($0.02/hr, 1 core , 613 MB RAM, 1-2  ECU)')
	print('                m1.large    ($0.34/hr, 2 cores, 7.5 GB RAM, 4    ECU)')
	print('                m1.xlarge   ($0.68/hr, 4 cores, 15  GB RAM, 8    ECU)')
	print('                c1.medium   ($0.17/hr, 2 cores, 1.7 GB RAM, 5    ECU)')
	print('                c1.xlarge   ($0.68/hr, 8 cores, 7   GB RAM, 20   ECU)')
	print('                cc1.4xlarge ($1.60/hr, 8 cores, 23  GB RAM, 33.5 ECU)')
	print('')
	print('   addSpotNodes [instance-type] [number] [maxprice] <group> -- adds one or')
	print('                more spot node request to the request queue. Apart from the')
	print('                maximum price parameter, the arguments are as above.')
	print('')
	print('   spotPrices [instance-type] -- get a 1-day history of spot instance')
	print('                prices for the specified type of machine (e.g. c1.xlarge)')
	print('')
	print('   terminate [nodename] -- Terminates a certain node specified by name. ')
	print('                Note that Amazon will still bill you for any partial ')
	print('                accrued hours!')
	print('')
	print('   terminateAll <group> -- Terminates *all* currently running nodes of the')
	print('                specified group. Note that Amazon will still bill you for ')
	print('                any partial accrued hours!')
	print('')
	print('   cancelSpot <spot-id> -- Cancels an active spot node request. The')
	print('                required ID can be obtained using the "status" command.')
	print('')
	print('   cancelAllSpot -- Cancels all active spot node requests.')
	print('')
	print('   login [nodename] -- Opens a shell on any running node, whose name')
	print('                can be obtained using the status command. All mitsuba-') 
	print('                related data (log files, etc) can be found in "/mnt".') 
	print('')
	print('   install <group> -- Fetch and install the most recent version of Mitsuba')
	print('                on all running compute nodes of the specied group.') 
	print('')
	print('   start <group> -- Start/restart the Mitsuba server on all currently running')
	print('                nodes of the specified group. The launcher creates a setup')
	print('                where all EC2 nodes are linked to each other on the fast ')
	print('                internal Amazon network, and any communication with the ')
	print('                outside world happens through a "head" node, whose address is')
	print('                printed at the end.')
	print('')
	print('   systemLoad <group> -- Prints the system load (1, 5 and 15-minute averages)')
	print('                for each machine in the specified group')
	print('')
	print('   runCommand "cmd args" <group> -- Runs the specified command on each machine')
	print('                in the specified group. Note that it has to be in quotes')
	print('')
	print('   syncData [prefix] <group> -- Downloads a file from the registered S3 bucket')
	print('               so that it is available to any rendering jobs. The download is')
	print('               simultaneously performed on all nodes of the specified group. ')
	print('               Note that only a prefix must be specified -- e.g. when there are')
	print('               20 data files starting with myData_*, running \"syncData myData\"')
	print('               will fetch them all. Any previously downloaded data will')
	print('               not be downloaded again. Note that you must ensure that a .s3cfg')
	print('               file is located in the same directory as this script, which ')
	print('               contains the S3 access credentials required by the "s3cmd"')
	print('               utility.')
	print('')
	print('  The usual order of these is: addNodes, install, start,')
	print('  and optionally syncData. At the end, dont\'t forget to run terminateAll.')
	print('  Many commands accept an optional <group> argument -- this can be used')
	print('  to set up multiple independent clusters on the same AWS account. If no')
	print('  group is specified, the default group ("default") is assumed.')
	print('')
	sys.exit(0)

conn = boto.ec2.connect_to_region(AWS_REGION,
	aws_access_key_id = AWS_ACCESS_KEY_ID, 
	aws_secret_access_key = AWS_SECRET_ACCESS_KEY)

if sys.argv[1] == 'addNodes':
	if len(sys.argv) == 5:
		addNodes(sys.argv[2], int(sys.argv[3]), sys.argv[4])
	elif len(sys.argv) == 4:
		addNodes(sys.argv[2], int(sys.argv[3]), 'default')
	else:
		print('addNodes: Invalid number of arguments!')
elif sys.argv[1] == 'addSpotNodes':
	if len(sys.argv) == 6:
		addSpotNodes(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
	elif len(sys.argv) == 5:
		addSpotNodes(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), 'default')
	else:
		print('addNodes: Invalid number of arguments!')
elif sys.argv[1] == 'status':
	if len(sys.argv) == 2:
		status()
	else:
		print('status: Invalid number of arguments!')
elif sys.argv[1] == 'terminate':
	if len(sys.argv) == 3:
		terminate(sys.argv[2])
	else:
		print('terminate: Invalid number of arguments!')
elif sys.argv[1] == 'terminateAll':
	if len(sys.argv) == 2:
		terminateAll('default')
	elif len(sys.argv) == 3:
		terminateAll(sys.argv[2])
	else:
		print('terminateAll: Invalid number of arguments!')
elif sys.argv[1] == 'cancelSpot':
	if len(sys.argv) == 3:
		cancelSpot(sys.argv[2])
	else:
		print('cancelSpot: Invalid number of arguments!')
elif sys.argv[1] == 'cancelAllSpot':
	if len(sys.argv) == 2:
		cancelAllSpot()
	else:
		print('cancelAllSpot: Invalid number of arguments!')
elif sys.argv[1] == 'install':
	if len(sys.argv) == 2:
		install('default')
	elif len(sys.argv) == 3:
		install(sys.argv[2])
	else:
		print('install: Invalid number of arguments!')
elif sys.argv[1] == 'syncData':
	if len(sys.argv) == 3:
		syncData(sys.argv[2], 'default')
	elif len(sys.argv) == 4:
		syncData(sys.argv[2], sys.argv[3])
	else:
		print('syncData: Invalid number of arguments!')
elif sys.argv[1] == 'start':
	if len(sys.argv) == 2:
		start('default')
	elif len(sys.argv) == 3:
		start(sys.argv[2])
	else:
		print('start: Invalid number of arguments!')
elif sys.argv[1] == 'login':
	if len(sys.argv) == 3:
		login(sys.argv[2])
	else:
		print('login: Invalid number of arguments!')
elif sys.argv[1] == 'spotPrices':
	if len(sys.argv) == 3:
		spotPrices(sys.argv[2])
	else:
		print('spotPrices: Invalid number of arguments!')
elif sys.argv[1] == 'systemLoad':
	if len(sys.argv) == 2:
		systemLoad('default')
	elif len(sys.argv) == 3:
		systemLoad(sys.argv[2])
	else:
		systemLoad('systemLoad: Invalid number of arguments!')
elif sys.argv[1] == 'runCommand':
	if len(sys.argv) == 3:
		runCommand(sys.argv[2], 'default')
	elif len(sys.argv) == 4:
		runCommand(sys.argv[2], sys.argv[3])
	else:
		print('runCommand: Invalid number of arguments!')
elif sys.argv[1] == 'regions':
	for r in conn.get_all_regions():
		print(r.name)
else:
	print('Unsupported command (run without parameters for an overview)!')
