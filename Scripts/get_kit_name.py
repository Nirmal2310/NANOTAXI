from minknow_api.manager import Manager

import re

client = list(Manager(host='localhost', port=9502).flow_cell_positions())[0].connect()

kit_name = str(client.protocol.get_run_info().user_info.kit_info).strip()

kit_name = re.sub(r'sequencing_kit:\s*"?(.*?)"?$', r'\1', kit_name)