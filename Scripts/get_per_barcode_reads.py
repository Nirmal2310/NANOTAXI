from minknow_api.manager import Manager

import minknow_api.statistics_pb2

import pandas as pd

# Connect to the manager

client = list(Manager(host='localhost', port=9502).flow_cell_positions())[0].connect()

# Get the Acquisition Run ID

run_id = client.acquisition.get_acquisition_info().run_id

data = []

data_stream = client.statistics.stream_acquisition_output(acquisition_run_id = client.acquisition.get_acquisition_info().run_id,
                                            data_selection=minknow_api.statistics_pb2.DataSelection(start=-1, end=-1),
                                            split=minknow_api.statistics_pb2.AcquisitionOutputSplit(barcode_name=True,
                                                                                                    lamp_barcode_id=True,
                                                                                                    lamp_target_id=True))

for filter_groups in data_stream:
    for filter_group in filter_groups.snapshots:
        for barcode in filter_group.filtering:
            for snapshot in filter_group.snapshots:
                data.append({
                    "Barcode": barcode.barcode_name,
                    "Counts": snapshot.yield_summary.basecalled_pass_read_count
                })
            
df = pd.DataFrame(data)