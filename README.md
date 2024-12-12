# AP_scripts_peterslab
Scripts written after opening lab

*NOTE*: all functions are in `+ap` namespace, so are called as `ap.function`

#### Loading recordings
`load_recording` - calls: 
- `load_timelite/bonsai/mousecam/widefield/ephys`, based on what's available and workspace structure `load_parts`
- parses wheel with `ap.parse_wheel` into `wheel velocity` and binary `wheel_move` vectors

#### Ephys
##### Preprocessing
- `preprocess_neuropixels(animal,day)` - saves metadata, runs PyKilosort, translates spike times into Open Ephys timestamps
##### Analysis
- `spike_cg` - auto/cross-correlation

#### Widefield 
`wf_align(im_unaligned,animal,day,align_type_master_align)`:
- create and save new day alignment: `ap.align_widefield([],animal,[],'new_days')` (uses all wf days, or specify `im_unaligned/days` for selection)
- create and save new animal alignment: `ap.align_widefield(day-aligned retinotopy,animal,[],'new_animal')` (aligns day-aligned average animal retinotopy to master retinotopy)

- `wf_roi`: grab fluorescence trace in ROI
- `wf_draw`: draw CCF borders/grid/point/areas over aligned widefield
- `wf_reflect`: flips L/R on aligned widefield image
- `wf_load_raw`: load/scroll raw frames
- `wf_corrviewer`: view seed pixel correlations in widefield


#### Data exploration
- `expscroll` - scroll through data

#### Plotting
- `ccf_draw` - draw CCF flat and 3D views, draw areas. plot probe trajectories
- `plot_probe_positions` - draw all probe positions (NTE and histology) for animal

#### Histology
- `where2slice` - show anterior/posterior extents of coronal slices from saved NTE trajectories


