function ise_json_to_hdf5(json_dir, hdf5_dir, start_year, end_year)
    time_stamps, datasets = ise_json_to_dict("public_power", "unix_seconds", "production_types", json_dir, start_year=start_year, end_year=end_year);
    
    time_stamps, datasets = ise_json_to_dict("total_power", "unix_seconds","production_types", json_dir, start_year=start_year, end_year=end_year);

    time_stamps, datasets = ise_json_to_dict("cbpf", "unix_seconds", "countries", json_dir, start_year=start_year, end_year=end_year);
    cbpf_hdf5_path = save_ise_data_to_hdf5(hdf5_dir,"cbpf", time_stamps, datasets, start_year, end_year)

    installed_power_to_hdf5(json_dir, hdf5_dir)
end
