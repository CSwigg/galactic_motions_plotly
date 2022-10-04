from webbrowser import get
import pandas as pd
import numpy as np
import movie_stars
import traces
import warnings
from astropy.io import fits

warnings.filterwarnings("ignore")
import plotly.graph_objects as go

import dash
from dash import html
from dash import dcc
from dash.dependencies import Input, Output, State

app = dash.Dash(__name__)
server = app.server
# ------------------------------------------------------------------------------
app.layout = html.Div([
    #html.Br(),
    dcc.Graph(id='map_2d',
              style={
                  'width': '82vh',
                  'height': '98vh',
                  'display': 'inline-block'
              }),
    dcc.Graph(id='map_3d',
              style={
                  'width': '90vh',
                  'height': '98vh',
                  'display': 'inline-block'
              })
])


def get_init_3d_fig_info():



    data_list, time_to_integrate = traces.traces()
    movie = movie_stars.Movie(movie_stars=data_list,
                                time=time_to_integrate,
                                movie_save_path=None,
                                center_star='Sun',
                                plot_dimensions=3,
                                annotations=None,
                                xyz_ranges=[4000, 4000, 500],
                                camera_follow=True,
                                center_galactic=False,
                                add_extra_traces=True)
    fig = movie.make_movie()
    return fig


def change_fig_symbols(fig, selectedData):
    point_names = []
    
    for point in selectedData['points']:
        point_names.append(point['hovertext'])

    new_frame_list = []
    for i, frame in enumerate(fig['frames']):
        data_list = []
        for j, data in enumerate(frame['data']):

            names = np.asarray(data['hovertext'])
            symbols = data['marker']['symbol']

            if (data['name'] != 'Zucker Clouds') & (data['name'] != 'Constant GC Lines'):

                if type(symbols) == str:
                    symbols = [symbols]

                symbols = np.asarray(symbols)
                symbols[np.where(np.isin(names, point_names))] = 'x'
                data['marker']['symbol'] = symbols
            data_list.append(data)
        frame['data'] = data_list
        new_frame_list.append(frame)
    fig['frames'] = new_frame_list
    return fig


@app.callback(Output(component_id='map_3d', component_property='figure'),
              Output(component_id='map_2d', component_property='figure'),
              Input(component_id='map_2d', component_property='selectedData'))
def callback_function(selectedData):

    fig = get_init_3d_fig_info()

    if selectedData is not None:
        if len(selectedData['points']) == 0:
            return dash.no_update, dash.no_update
        fig = change_fig_symbols(fig=fig, selectedData=selectedData)

    data_list = fig['frames'][0]['data']
    new_data_list = []
    for i, data in enumerate(data_list):
        if 'hovertext' in data:
            hovertext = data['hovertext']
        else:
            hovertext = None
        marker = data['marker']

        scatter = go.Scatter(x=data['x'],
                             y=data['y'],
                             mode='markers',
                             marker=dict(color=marker['color'],
                                         symbol=marker['symbol'],
                                         opacity=marker['opacity'],
                                         size=marker['size'],
                                         line=dict(width=0.)),
                             line=dict(width=0.00001),
                             name=data['name'],
                             hovertext=hovertext,
                             customdata=hovertext)
        new_data_list.append(scatter)
    fig_2d = {}
    fig_2d['data'] = new_data_list
    fig_2d['layout'] = go.Layout(template = 'plotly_dark', 
                                xaxis=dict(autorange=False,
                                            range=[-4000, 4000],
                                            showgrid=False,
                                            zeroline=False,
                                            nticks=3,
                                            linecolor='white',
                                            linewidth=3,
                                            title=dict(text='X (pc)')),
                                 yaxis=dict(autorange=False,
                                            range=[-4000, 4000],
                                            showgrid=False,
                                            zeroline=False,
                                            nticks=3,
                                            linecolor='white',
                                            linewidth=3,
                                            title=dict(text='Y (pc)')))
    fig_2d['layout']['legend'] = dict(x=0,
                                        y=1.,
                                        font=dict(size=13,
                                                    family='Courier New'),
                                        itemsizing='constant',
                                        bgcolor='rgba(0,0,0,0)')
    return go.Figure(fig), go.Figure(fig_2d)


if __name__ == '__main__':
    app.run_server(debug=True)
