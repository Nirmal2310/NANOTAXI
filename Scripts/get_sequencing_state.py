from minknow_api.manager import Manager

# Connect to the manager

client = list(Manager(host='localhost', port=9502).flow_cell_positions())[0].connect()

while True:
    try:
        # Get the sequencing state
        phase = client.protocol.get_current_protocol_run().phase
        if(phase == 1):
            status = print("Initializing")
        elif(phase == 2):
            status = print("Sequencing")
        elif(phase == 4):
            status = print("Pore Scanning")
    except Exception as e:
        status = print("Completed")
        break


