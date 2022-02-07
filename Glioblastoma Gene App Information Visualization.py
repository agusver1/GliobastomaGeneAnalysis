import json
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import dash_cytoscape as cyto #pip install dash-cytoscape
import networkx as nx
import numpy as np
import matplotlib.cm as cm 
import dash_daq as daq #pip install dash_daq
import plotly.express as px



app = dash.Dash(__name__)
server = app.server

#User .txt path load:
    
nodeload='C:/Users/agusv/Desktop/Estudio/LDATA2010-Information Visualization/PROYECTO/BIOGRID-PROJECT-glioblastoma_project-GENES.projectindex.txt'
edgeload='C:/Users/agusv/Desktop/Estudio/LDATA2010-Information Visualization/PROYECTO/BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS.tab3.txt'
chemicalload='C:/Users/agusv/Desktop/Estudio/LDATA2010-Information Visualization/PROYECTO/BIOGRID-PROJECT-glioblastoma_project-CHEMICALS.chemtab.txt'
ptmload='C:/Users/agusv/Desktop/Estudio/LDATA2010-Information Visualization/PROYECTO/BIOGRID-PROJECT-glioblastoma_project-PTM.ptmtab.txt'



#Load Glioblastoma Node .txt file
read_nodes=pd.read_csv(nodeload,sep='\t')
datanames=read_nodes.columns[1:]
read_nodes = read_nodes.set_index('#BIOGRID ID')
read_nodes = read_nodes.T.to_dict('list')
# ###################### DATA PREPROCESSING ######################
#Load Glioblastoma Edges .txt file
read_edges=pd.read_csv(edgeload,sep='\t')
a=list(read_edges['BioGRID ID Interactor A'][1:])
b=list(read_edges['BioGRID ID Interactor B'][1:])
edges=list(zip(a,b))

#Remove duplicate connections:
edges=set(edges)

#Load Chemical Data:
chem_data=pd.read_csv(chemicalload,sep='\t')
chemnames=chem_data.columns[1:]
chem_data = chem_data.set_index('BioGRID Gene ID')
chem_data = chem_data.T.to_dict('list')

#Load PTM Data:
ptm_data=pd.read_csv(ptmload,sep='\t')
ptmnames=ptm_data.columns[1:]
ptm_data = ptm_data.set_index('BioGRID ID')
ptm_data = ptm_data.T.to_dict('list')

#To start computing the algorithms. We use NetworkX.

G = nx.Graph()
G.add_edges_from(edges)

#Function to normalize the values obtaines from the NetworkX algorithms.
def algorithm_normalize(weights):
    values=[]
    for node in weights:
        values.append(weights[node])
    listvalues=[(g-min(values))/(max(values)-min(values)) for g in values]
    return listvalues

nodict={}

#Load Extra Cytoscape Layouts:
cyto.load_extra_layouts()

#Turning NetworkX algorithm values to Hexagesimal to color the nodes:
def cmap_at_scale(cmap, bins, nodes):
    # cmap:  mpl or cmocean colormap 
    # bins :  ascending ordered list  of floats in [0,1] that defines the scale for a plotly colorscale derived from cmap
    # returns   a list of hexcolors
    
    scale  = bins
    colors = cmap(scale)[:, :3]
    colors = (255*colors+0.5).astype(np.uint8)
    hexcolors = [f'#{c[0]:02x}{c[1]:02x}{c[2]:02x}' for c in colors]
    for i in range(len(hexcolors)):
        nodict[nodes[i]]=hexcolors[i]
    return nodict

#Charge nodes and edges:

nodes = []

tix=""
cy_edges = []
cy_nodes = []


#Counter of the connections each node has:

counter={}

#Add nodes and edges to its respective variables. Also add the information each one has in the DataSet.

for network_edge in edges:
    source=network_edge[0]
    target=network_edge[1]
    
    if source not in nodes:
        nodes.append(source)
        counter[str(source)]=1
        cy_nodes.append({"data": {"id": source, "label": "ID #" + str(source),'size':1}})

    else:
        counter[str(source)]+=1
        num=list(nodes).index(source)
        cy_nodes[num]['data']['size']=int(counter[str(source)])
    
    if target not in nodes:
        nodes.append(target)
        counter[str(target)]=1
        cy_nodes.append({"data": {"id": target, "label": "ID #" + str(target),'size':1}})

    else:
        counter[str(target)]+=1      
        num1=list(nodes).index(target)
        cy_nodes[num1]['data']['size']=int(counter[str(target)])
        
#Create edge list and prevent self-connections or multiple connections between the same nodes:
    if source != target:
        cy_edges.append({
                            'data': {
                                'source': source,
                                'target': target,
                            }
                        })

#Node self information charge in dictionary:
tix={}
for node in nodes:
    tix[str(node)]=""
    if node in read_nodes.keys():
        for i in range(len(datanames)):
            tix[str(node)]= tix[str(node)]+ str(datanames[i])+':' + '\t' + str(read_nodes[node][i])+ "\n"
    else:
        tix[str(node)]="\n"+"No more information available about this gene in the database."

#Chemical information for treatment:

cheminfo={}
for node in nodes:
    cheminfo[str(node)]=""
    if node in chem_data.keys():
        for i in range(len(chemnames)):
            cheminfo[str(node)]= cheminfo[str(node)]+ str(chemnames[i])+':' + '\t' + str(chem_data[node][i])+ "\n"
    else:
        cheminfo[str(node)]="\n"+"No more information available about this gene in the Chemicals database."

#PTM information for treatment:

ptminfo={}
for node in nodes:
    ptminfo[str(node)]=""
    if node in ptm_data.keys():
        for i in range(len(ptmnames)):
            ptminfo[str(node)]= ptminfo[str(node)]+ str(ptmnames[i])+':' + '\t' + str(ptm_data[node][i])+ "\n"
    else:
        ptminfo[str(node)]="\n"+"No more information available about this gene in the Chemicals database."
#Stylesheet dash is gonna use.
            
default_stylesheet = [
    {
        "selector": 'node',
        'style': {
            "opacity": 1,
        }
    },
    {
        "selector": 'edge',
        'style': {
            "curve-style": "bezier",
            "opacity": 1
        }
    },
]

styles = {
    'json-output': {
        'overflow-y': 'scroll',
        'height': 'calc(50% - 25px)',
        'border': 'thin lightgrey solid'
    },
    'tab': {
        'height': 'calc(98vh - 105px)'
    }
}

#We start the Dash App.

app.layout = html.Div([
    html.Div(className='eight columns', children=[
        cyto.Cytoscape(
            id='cytoscape',
            elements=[],
            style={
                'height': '95vh',
                'width': '100%'
            }
        )
    ]),

    html.Div(className='four columns', children=[
        dcc.Tabs(id='tabs', children=[
            dcc.Tab(label='Control Panel', children=[
                html.Div([html.H2('Layout chosen:'),
                dcc.Dropdown(
                    id='dropdown-layout',
                    options=[
                        {'label': 'Random', 'value': 'random'},
                        {'label': 'Grid', 'value': 'grid'},
                        {'label': 'Circle', 'value': 'circle'},
                        {'label': 'Concentric', 'value': 'concentric'},
                        {'label': 'Breadthfirst', 'value': 'breadthfirst'},
                        {'label': 'Cose', 'value': 'cose'},
                        {'label': 'Cose Bilkent', 'value': 'cose-bilkent'}
                        
                    ],
                    value='grid',
                    clearable=False
                )]),
                html.Div([html.H2('Color ranking by Network Algorithms:'),
                dcc.Dropdown(
                    id='dropdown-algorithms',
                    options=[
                        {'label': 'Degree Centrality', 'value': 'nx.degree_centrality(G)'},
                        {'label': 'Closeness Centrality', 'value': 'nx.closeness_centrality(G)'},
                        {'label': 'Betweenness centrality', 'value': 'nx.betweenness_centrality(G)'},
                        {'label': 'Clustering Coefficient', 'value': 'nx.clustering(G)'},
                        {'label': 'Google PageRank', 'value': 'nx.pagerank(G)'}
                    ],
                    value='nx.degree_centrality(G)',
                    clearable=False
                )]),
                html.Div([html.H2('Node shape selection:'),
                dcc.Dropdown(
                    id='dropdown-node-shape',
                    value='ellipse',
                    clearable=False,
                    options=[
                        {'label': 'Ellipse', 'value': 'ellipse'},
                        {'label': 'Triangle', 'value': 'triangle'},
                        {'label': 'Rectangle', 'value': 'rectangle'},
                        {'label': 'Diamond', 'value': 'diamond'},
                        {'label': 'Pentagon', 'value': 'pentagon'},
                        {'label': 'Hexagon', 'value': 'hexagon'},
                        {'label': 'Heptagon', 'value': 'heptagon'},
                        {'label': 'Octagon', 'value': 'octagon'},
                        {'label': 'Star', 'value': 'star'},
                        {'label': 'Polygon', 'value': 'polygon'}
                    ]
                )]),
                html.Div([html.H2('Color map selection:'),
                dcc.Dropdown(
                    id='dropdown-colors',
                    value='cm.Wistia',
                    clearable=False,
                    options=[
                        {'label': 'Wistia', 'value': 'cm.Wistia'},
                        {'label': 'Jet', 'value': 'cm.jet'},
                        {'label': 'Binary', 'value': 'cm.binary'},
                        {'label': 'Greens', 'value': 'cm.Greens'}
                    ]
                )]),
                html.Div([html.H2('Number of connections per node more than:'),
                daq.NumericInput(
                    id='my-numeric-input-1',
                    value=150
                ),
                html.Div(id='numeric-input-output-1')]),
                html.Div([html.H2('Follower and Following node colors:'),
                dcc.Input(
                    id='input-follower-color',
                    type='text',
                    value='#0074D9',
                ),

                dcc.Input(
                    id='input-following-color',
                    type='text',
                    value='#FF4136',
                )]),
            ]),
            dcc.Tab(label='Information', children=[
                html.Div(style=styles['tab'], children=[

                    html.P('Node selected:',style={
                              "fontSize": 18,
                              "font-family": 'georgia',
                              "font-weight": "bold",
                              "text-align": "left",
                            }),
                    html.Pre(
                        id='tap-node-information',
                        style={
                                "fontSize": 15,
                                "font-family": 'georgia',
                                "text-align": "left",
                            }),
                    
                    
                    html.P('Node targets:',style={
                              "fontSize": 18,
                              "font-family": 'georgia',
                              "font-weight": "bold",
                              "text-align": "left",
                            }),
                    html.Pre(
                        id='tap-node-information-target',
                        style={
                                "fontSize": 15,
                                "font-family": 'georgia',
                                "text-align": "left",
                            }
                    ), 
                    html.P('Node Chemical Information:',style={
                              "fontSize": 18,
                              "font-family": 'georgia',
                              "font-weight": "bold",
                              "text-align": "left",
                            }),
                    html.Pre(
                        id='tap-node-information-chemical',
                        style={
                                "fontSize": 15,
                                "font-family": 'georgia',
                                "text-align": "left",
                            }
                    ),
                    html.P('Post-Translational Node Modifications:',style={
                              "fontSize": 18,
                              "font-family": 'georgia',
                              "font-weight": "bold",
                              "text-align": "left",
                            }),
                    html.Pre(
                        id='tap-node-information-ptm',
                        style={
                                "fontSize": 15,
                                "font-family": 'georgia',
                                "text-align": "left",
                            }
                    ),
                ])
            ]),   
            dcc.Tab(label='JSON', children=[
                html.Div(style=styles['tab'], children=[
                    html.P('Node Object JSON:'),
                    html.Pre(
                        id='tap-node-json-output',
                        style=styles['json-output']
                    ),
                    html.P('Edge Object JSON:'),
                    html.Pre(
                        id='tap-edge-json-output',
                        style=styles['json-output']
                    )
                ])
            ]),
            dcc.Tab(label='Histogram', children=[
                html.Div(style=styles['tab'], children=[
                    html.P('Graph analysis with the chosen method:'),
                    dcc.Graph(id="histogram"),
                    dcc.Graph(id="heatmap")

                ])
            ]),
        ]),
    ])
])


@app.callback(Output('tap-node-json-output', 'children'),
              [Input('cytoscape', 'tapNode')])
def display_tap_node(data):
    return json.dumps(data, indent=2)

#Information regarding the node sources in the DataBase.

@app.callback(Output('tap-node-information', 'children'),
              [Input('cytoscape', 'tapNode')])
def displayTapNodeData(node):
    
    if not node:
        text=''
        return text
    else:
        ashe=node['data']['id']
        texto="ID:"+'\t'+str(ashe)+"\n"+'NUMBER OF INTERACTIONS:'+'\t'+ str(node['data']['size']) + "\n"
        texto=texto + tix[ashe]
        text = [html.P(texto)]
    return text

#Information regarding the targets in the nodes database.

@app.callback(Output('tap-node-information-target', 'children'),
              [Input('cytoscape', 'tapNode')])
def displayTapTarget(node):
    if not node:
        text=''
        return text
    elif node['edgesData']==[]:
        text='No connections detected for this node'
        return text
    else:
        for i in range(len(node['edgesData'])):
            ashe=node['edgesData']
            ashei=ashe[i]['target']
            texto="ID:"+'\t'+str(ashei)+"\n"+'NUMBER OF INTERACTIONS:'+'\t'+ str(counter[ashei]) + "\n"
            texto=texto + tix[ashei]
            text = [html.P(texto)]
    return text

#Information regarding the Chemical reactions.

@app.callback(Output('tap-node-information-chemical', 'children'),
              [Input('cytoscape', 'tapNode')])
def displayTapChemical(node):
    if not node:
        text=''
        return text
    else:
        ashe=node['data']['id']
        texto="This node contains chemical information in the Database" + "\n"
        texto=texto + cheminfo[ashe]
        text = [html.P(texto)]
    return text

#Information regarding the PTM.

@app.callback(Output('tap-node-information-ptm', 'children'),
              [Input('cytoscape', 'tapNode')])
def displayTapPtm(node):
    if not node:
        text=''
        return text
    else:
        ashe=node['data']['id']
        texto="This node contains chemical information in the Database" + "\n"
        texto=texto + ptminfo[ashe]
        text = [html.P(texto)]
    return text

#Show JSON edge information in the tab.

@app.callback(Output('tap-edge-json-output', 'children'),
              [Input('cytoscape', 'tapEdge')])
def display_tap_edge(data):
    return json.dumps(data, indent=2)

#Update nodes and edges as the user applies connection filters.

@app.callback(Output('cytoscape', 'elements'),
              [Input('my-numeric-input-1', 'value'),
               Input('dropdown-algorithms', 'value'),
               Input('dropdown-colors', 'value')])

def update_cytoscape_elements(value,algorithm,color):
    eliminator=[]
    nodenum=len(cy_nodes)
    cy_edges2=[]
    measure=eval(algorithm)
    filters=algorithm_normalize(measure)
    algorithm_weight=cmap_at_scale(eval(color), filters, nodes)
    
    for node in nodes:
        num=list(nodes).index(node)
        cy_nodes[num]['data']['algorithm']=algorithm_weight[node]
    
    for i in range(nodenum):
        if cy_nodes[i]['data']['size'] < value:
            eliminator.append(cy_nodes[i]['data']['id'])
    cy_nodes2 = [x for x in cy_nodes if x['data']['id'] not in eliminator]
    
    for x in cy_edges:
        if x['data']['source'] not in eliminator:
            if x['data']['target'] not in eliminator:
                cy_edges2.append(x)
    return cy_edges2 + cy_nodes2

#Update Layout as the user chooses a value in the dropdown.

@app.callback(Output('cytoscape', 'layout'),
              [Input('dropdown-layout', 'value')])
def update_cytoscape_layout(layout):
        return {'name': layout}

#Change the stylesheet as the user taps a node.

@app.callback(Output('cytoscape', 'stylesheet'),
              [Input('cytoscape', 'tapNode'),
               Input('input-follower-color', 'value'),
               Input('input-following-color', 'value'),
               Input('dropdown-node-shape', 'value')])


def generate_stylesheet(node, follower_color, following_color, node_shape):

    if not node:
        default_stylesheet = [
                {
                    "selector": 'node',
                    'style': {
                        'background-color': 'data(algorithm)',
                    }
                },
                {
                    "selector": 'edge',
                    'style': {
                        "curve-style": "bezier",
                        "opacity": 1
                    }
                }]
        return default_stylesheet
        
    stylesheet = [{
        "selector": 'node',
        'style': {
            'opacity': 0.3,
            'shape': node_shape,
        }
    }, {
        'selector': 'edge',
        'style': {
            'opacity': 0.2,
            "curve-style": "bezier",
        }
    }, {
        "selector": 'node[id = "{}"]'.format(node['data']['id']),
        "style": {
            'background-color': '#B10DC9',
            "border-color": "purple",
            "border-width": 2,
            "border-opacity": 1,
            "opacity": 1,

            "label": "data(label)",
            "color": "#B10DC9",
            "text-opacity": 1,
            "font-size": 12,
            'z-index': 9999
        }
    }]

    for edge in node['edgesData']:
        if edge['source'] == node['data']['id']:
            stylesheet.append({
                "selector": 'node[id = "{}"]'.format(edge['target']),
                "style": {
                    'background-color': following_color,
                    'opacity': 0.9
                }
            })
            stylesheet.append({
                "selector": 'edge[id= "{}"]'.format(edge['id']),
                "style": {
                    "mid-target-arrow-color": following_color,
                    "mid-target-arrow-shape": "vee",
                    "line-color": following_color,
                    'opacity': 0.9,
                    'z-index': 5000
                }
            })

        if edge['target'] == node['data']['id']:
            stylesheet.append({
                "selector": 'node[id = "{}"]'.format(edge['source']),
                "style": {
                    'background-color': follower_color,
                    'opacity': 0.9,
                    'z-index': 9999
                }
            })
            stylesheet.append({
                "selector": 'edge[id= "{}"]'.format(edge['id']),
                "style": {
                    "mid-target-arrow-color": follower_color,
                    "mid-target-arrow-shape": "vee",
                    "line-color": follower_color,
                    'opacity': 1,
                    'z-index': 5000
                }
            })

    return stylesheet

@app.callback([Output("histogram", "figure"),
               Output("heatmap",'figure')],
              [Input('dropdown-algorithms', 'value')])

def display_color(algorithm):
    measure=eval(algorithm)
    filters=algorithm_normalize(measure)
    fig = px.histogram(filters, range_x=[0, 1], range_y=[0, 100],
                       labels={
                           'value':'Node values according to algorithm',
                           'count':'Frequency',
                           'variable':'Node values'
                           }, title="Histogram due to Algorithm Values")
    degrees = [len(list(G.neighbors(n))) for n in G.nodes()]
    data_tuples = list(zip(degrees,filters))
    data=pd.DataFrame(data_tuples, columns=['Degrees','Value'])
    fig2= px.scatter(x=data['Degrees'],y=data['Value'],color=data['Value'],
                     labels={
                           'x':'Number of Neighbours',
                           'y':'Value in Algorithm'},
                     title="Number of Neighbours vs Value in Algorithm")

    return fig, fig2



if __name__ == '__main__':
    app.run_server(debug=True)
