from minknow_api.manager import Manager

# Connect to the Sequencing Run

client = list(Manager(host='localhost', port=9502).flow_cell_positions())[0].connect()

# Get Pore Scanning Output

mux_scan_data = client.acquisition.get_acquisition_info().bream_info.mux_scan_results

# Print the number of active pores

print(mux_scan_data[0].counts.get("single_pore"))