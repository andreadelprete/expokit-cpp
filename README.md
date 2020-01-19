# expokit-cpp

Trying to make it compile with f2c requirement

1.  Install f2c via apt
2.  Add the followings to the env var `PKG_CONFIG_PATH`, according to your system, but stuff should match:
    1. `/opt/openrobots/share/pkgconfig`
    2. `/opt/openrobots/lib/pkgconfig`
    3. `<any path sensible to you>`
3.  Add the `f2c.pc` file to `<any path sensible to you>`
