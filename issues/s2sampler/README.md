see https://forum.step.esa.int/t/large-time-difference-when-using-core-or-s2tbx-resamplers-for-s2-data-with-view-angles-from/11791/3

run example:
```bash
/usr/bin/time -f "RSS=%M elapsed=%E cpu.sys=%S .user=%U" python3 test2_snappy_resampling.py 'core' S2A_MSIL1C_20170210T082051_N0204_R121_T33HYD_20170210T083752.SAFE
```