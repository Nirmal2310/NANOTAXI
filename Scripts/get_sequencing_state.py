from minknow_api.manager import Manager

try:

	# Connect to the manager
	client = list(Manager(host='localhost', port=9502).flow_cell_positions())[0].connect()
	# Get the sequencing state
	phase = client.protocol.get_current_protocol_run().phase

	if(phase == 1):
		status = "Initializing"
	elif(phase == 2):
		status = "Sequencing"
	elif(phase == 4):
		status = "Pore Scanning"

except Exception as e:
	print("No Sequencing Run is Going On.")

#logging.get_absl_handler().python_handler.stream = sys.stdout

