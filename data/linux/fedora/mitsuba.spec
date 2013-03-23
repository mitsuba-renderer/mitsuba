Name:		mitsuba
Version:	0.4.4
Release:	1%{?dist}
Summary:	Mitsuba renderer
Group:		Applications/Graphics
License:	GPL-3
URL:		http://www.mitsuba-renderer.org
Source:		%{name}-%{version}.tar.gz
BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)
BuildRequires:	gcc-c++ scons boost-devel qt4-devel OpenEXR-devel xerces-c-devel python-devel glew-devel collada-dom-devel eigen3-devel
Requires:	boost qt4 OpenEXR-libs xerces-c python libGLEWmx collada-dom
%description
Mitsuba is an extensible rendering framework written in portable C++. It implements unbiased as well as biased techniques and contains heavy optimizations targeted towards current CPU architectures.

The program currently runs on Linux, MacOS X and Microsoft Windows and makes use of SSE2 optimizations on x86 and x86_64 platforms. So far, its main use has been as a testbed for algorithm development in computer graphics, but there are many other interesting applications.

Mitsuba comes with a command-line interface as well as a graphical frontend to interactively explore scenes. While navigating, a rough preview is shown that becomes increasingly accurate as soon as all movements are stopped. Once a viewpoint has been chosen, a wide range of rendering techniques can be used to generate images, and their parameters can be tuned from within the program.
%package devel
Summary:	Mitsuba development files
Requires:	boost-devel qt4-devel OpenEXR-devel xerces-c-devel python-devel glew-devel collada-dom-devel
%description devel
This package contains the development headers, which are needed when
building custom plugins and other extensions for Mitsuba.
%prep
%setup -q
%build
cat build/config-linux-gcc.py | sed -e "s/\(boost_[^']*\)/\1-mt/g" > config.py
scons
%install
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT%{_libdir}
mkdir -p $RPM_BUILD_ROOT%{_bindir}
mkdir -p $RPM_BUILD_ROOT%{_libdir}/python2.7/lib-dynload
mkdir -p $RPM_BUILD_ROOT/usr/share/mitsuba/plugins
mkdir -p $RPM_BUILD_ROOT/usr/share/pixmaps
mkdir -p $RPM_BUILD_ROOT/usr/share/applications
mkdir -p $RPM_BUILD_ROOT/usr/include
strip dist/lib* dist/mtsgui dist/mitsuba dist/mtssrv dist/mtsutil dist/mtsimport
strip dist/plugins/* dist/python/*/*
cp dist/libmitsuba-*.so $RPM_BUILD_ROOT%{_libdir}
cp dist/mtsgui $RPM_BUILD_ROOT%{_bindir}
cp dist/mitsuba $RPM_BUILD_ROOT%{_bindir}
cp dist/mtssrv $RPM_BUILD_ROOT%{_bindir}
cp dist/mtsutil $RPM_BUILD_ROOT%{_bindir}
cp dist/mtsimport $RPM_BUILD_ROOT%{_bindir}
cp dist/python/2.7/mitsuba.so $RPM_BUILD_ROOT%{_libdir}/python2.7/lib-dynload
cp dist/plugins/* $RPM_BUILD_ROOT/usr/share/mitsuba/plugins
cp -Rdp dist/data $RPM_BUILD_ROOT/usr/share/mitsuba/data
cp src/mtsgui/resources/mitsuba48.png $RPM_BUILD_ROOT/usr/share/pixmaps
cp data/linux/mitsuba.desktop $RPM_BUILD_ROOT/usr/share/applications
cp -Rdp include/mitsuba $RPM_BUILD_ROOT/usr/include/mitsuba
%clean
rm -rf $RPM_BUILD_ROOT
%files
%defattr(-,root,root,-)
%{_libdir}/libmitsuba-*.so
%{_libdir}/python2.7/lib-dynload/mitsuba.so
%{_bindir}/*
/usr/share/pixmaps/mitsuba48.png
/usr/share/applications/mitsuba.desktop
/usr/share/mitsuba/*
%files devel
/usr/include/*
%changelog

* Thu Feb 28 2013 Wenzel Jakob <wenzel@cs.cornell.edu> 0.4.4%{?dist}
- Upgrade to version 0.4.4

* Tue Jan 29 2013 Wenzel Jakob <wenzel@cs.cornell.edu> 0.4.3%{?dist}
- Upgrade to version 0.4.3

* Wed Oct 31 2012 Wenzel Jakob <wenzel@cs.cornell.edu> 0.4.2%{?dist}
- Upgrade to version 0.4.2

* Wed Oct 10 2012 Wenzel Jakob <wenzel@cs.cornell.edu> 0.4.1%{?dist}
- Upgrade to version 0.4.1

* Thu Sep 27 2012 Wenzel Jakob <wenzel@cs.cornell.edu> 0.4.0%{?dist}
- Upgrade to version 0.4.0

* Fri Apr 13 2012 Wenzel Jakob <wenzel@cs.cornell.edu> 0.3.1%{?dist}
- Upgrade to version 0.3.1

* Mon Aug 15 2011 Wenzel Jakob <wenzel@cs.cornell.edu> 0.3.0%{?dist}
- First official Fedora Core build
