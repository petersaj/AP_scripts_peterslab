# AP_scripts_peterslab
Scripts written after opening lab

#### Loading recordings
`ap.load_recording` - calls: 
- `ap.load_timelite/bonsai/mousecam/widefield/ephys`, based on what's available and workspace structure `load_parts`
- parses wheel with `ap.parse_wheel` into `wheel velocity` and binary `wheel_move` vectors

#### Ephys
##### Preprocessing
`ap.preprocess_neuropixels(animal,day)` - saves metadata, runs PyKilosort, translates spike times into Open Ephys timestamps
##### Analysis
`ap.spike_cg` - auto/cross-correlation

#### Widefield 
`ap.wf_align(im_unaligned,animal,day,align_type_master_align)`:
- create and save new day alignment: `ap.align_widefield([],animal,[],'new_days')` (uses all wf days, or specify `im_unaligned/days` for selection)
- create and save new animal alignment: `ap.align_widefield(day-aligned retinotopy,animal,[],'new_animal')` (aligns day-aligned average animal retinotopy to master retinotopy)

`ap.wf_roi`: grab fluorescence trace in ROI
`ap.wf_draw`: draw CCF borders/grid/point over aligned widefield
`ap.wf_reflect`: flips L/R on aligned widefield image
`ap.wf_load_raw`: load/scroll raw frames
`ap.wf_corrviewer`: view seed pixel correlations in widefield


#### Data exploration
`ap.expscroll` - scroll through data

#### Plotting
`ccf_draw` - draw CCF flat and 3D views, draw areas. plot probe trajectories
`plot_probe_positions` - draw all probe positions (NTE and histology) for animal
