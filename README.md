# AP_scripts_peterslab
Scripts written after opening lab

#### Loading recordings
`ap.load_recording` - calls: 
- `ap.load_timelite/bonsai/mousecam/widefield/ephys`, based on what's available and workspace structure `load_parts`
- parses wheel with `ap.parse_wheel` into `wheel velocity` and binary `wheel_move` vectors

#### Ephys
##### Preprocessing
`ap.preprocess_neuropixels(animal,day)` - saves metadata, runs PyKilosort, translates spike times into Open Ephys timestamps

#### Widefield 
`ap.align_widefield(im_unaligned,animal,day,align_type_master_align)`:
- create and save new day alignment: `ap.align_widefield([],animal,[],'new_days')` (uses all wf days, or specify `im_unaligned/days` for selection)
- create and save new animal alignment: `ap.align_widefield(day-aligned retinotopy,animal,[],'new_animal')` (aligns day-aligned average animal retinotopy to master retinotopy)

`ap.draw_wf_ccf`: draw CCF borders over aligned widefield
`ap.reflect_widefield`: flips L/R on aligned widefield image

#### Data exploration
`ap.expscroll(wf_U,wf_V,wf_times,mousecam_fn,mousecam_times)` - scroll through widefield and/or mousecam
