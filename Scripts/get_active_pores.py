from minknow_api.manager import Manager

import pandas as pd

client = list(Manager(host='localhost', port=9502).flow_cell_positions())[0].connect()

channels = client.device.get_flow_cell_info().channel_count

if(channels==126):
    data = client.data.get_channel_states(first_channel = 1, last_channel = 126)
else:
    data = client.data.get_channel_states(first_channel = 1, last_channel = 2048)

channels_dict = []

for i in data:
   for channel in i.channel_states:
      channels_dict.append({
                    "Channel": channel.channel,
                    "Type": channel.state_name
                  })
   break

channels_df = pd.DataFrame(channels_dict)

channels_stat = channels_df['Type'].value_counts().reset_index()

channels_stat.columns = ['Type', 'Count']

pore_counts = channels_stat.loc[channels_stat['Type'].isin(['strand', 'adapter', 'pore']), 'Count'].sum()