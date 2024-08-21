from minknow_api.manager import Manager

# Connect to the manager 

client = list(Manager(host='localhost', port=9502).flow_cell_positions())[0].connect()

# Get the output directory for the current run

data_directory = client.protocol.get_run_info().output_path

print(data_directory)