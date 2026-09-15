%%%% CK ephys data

%% Plot recording locations

animals = {'HA016','HA017','HA018','HA019','HA020','AP036','AP037'};

animal_col = vertcat(ap.colormap('tube'),lines(7));
ccf_draw = ap.ccf_draw;
ccf_draw.draw_name('Caudoputamen');

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    probe_color = animal_col(curr_animal,:);
    ccf_draw.draw_probes_histology(animal,probe_color);
    drawnow;
end

