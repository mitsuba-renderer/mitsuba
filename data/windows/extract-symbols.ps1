Remove-Item -Recurse -Force symbols
mkdir symbols
Get-ChildItem build -include *.pdb -recurse | foreach ($_) { python dependencies\bin\symbolstore.py dependencies\bin\dump_syms_x86 symbols $_.fullname }