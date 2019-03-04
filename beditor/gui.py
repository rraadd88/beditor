import pandas as pd
import numpy as np
from os.path import dirname, basename, exists, abspath
from os import makedirs 
import PySimpleGUI as sg
import beditor
import multiprocessing
from beditor.lib.global_vars import multint2reg,nt2complement
from beditor.lib.io_strs import isstrallowed
from beditor.lib.io_dfs import *
from beditor.lib.global_vars import aminoacids
import subprocess
import logging

# vars beditor
nts=''.join(list(nt2complement.keys()))
# beditor io

def dropna(x):
    x_=[]
    for i in x:
        if not pd.isnull(i):
            x_.append(i)
    return x_
def unique(l,drop=None):
    l=[str(s) for s in l]
    if drop is not None:
        l=[s for s in l if s!=drop]
    return tuple(np.unique(l))
def unique_dropna(l): return dropna(unique(l,drop='nan'))
def get_dbepams():
    dbepams=pd.read_table(f'{dirname(abspath(__file__))}/data/dbepams.tsv',keep_default_na=False)
    dbepams['BE type']=dbepams.apply(lambda x : f"{x['nucleotide']}-{x['nucleotide mutation']}", axis=1)
    dbepams['editing window']=dbepams.apply(lambda x : f"{x['window start']}-{x['window end']}bp", axis=1)
    dbepams['BE type and PAM']=dbepams.apply(lambda x : ' '.join([f"{k}:{str(x[k])}" if k!='BE type' else f"{str(x[k])}" for k in ['BE type','PAM']]), axis=1)
    dbepams['BE name and editing window']=dbepams.apply(lambda x : ' '.join([f"{k}:{str(x[k])}" if k!='BE name' else f"{str(x[k])}" for k in ['method','editing window']]), axis=1)
    return dbepams.reset_index()

import yaml
def normalisestr(s):
    import re
    return re.sub('\W+','', s.lower()).replace('_','')

def loadcfginwinvals(win,vals): 
    for k in vals:
        if isinstance(vals[k],list): 
            for i in vals[k]:
                try:
                    win.FindElement(f"{k} {i}").Update(i)
                except:
                    logging.warning(f"{k} {i} element is not in window")
        else:
            win.FindElement(k).Update(vals[k])
    return win 

def get_mutation(vals2):
    mutation=[]
    for s in ['from','to']:
        keys=np.array([k for k in vals2 if str(k).startswith(s) and len(str(k))==len(str(s))+1])
        buls=np.array([vals2[k] for k in keys])
        mutation.append(keys[buls][0].replace(s,''))
    return '-'.join(mutation)


## coonvert vals to cfg
import numpy as np 
def guival2cfg(val,vals2):
#     val['dsubmap_preferred_path']=None if val['dsubmap_preferred_path']=='' else val['dsubmap_preferred_path']
#     val['mimetism_level']=None if val['mimetism_level']=='' else val['mimetism_level']

    cfg={}
    if not val["calculate beditor scores"]:
        cfg['step2ignore']=4
    else:
        cfg['step2ignore']=None        
    cfg['host']=val['Species name (Ensembl assembly)'].split(' (')[0]
    cfg['genomeassembly']=val['Species name (Ensembl assembly)'].split(' (')[1].replace(')','')
#     cfg['genomerelease']=int(val['genomerelease'].replace('genome release=',''))
    cfg['genomerelease']=95
    # cfg['BE and PAM']=[f"{val['BE name and editing window']} PAM:{pam}"]
    cfg['BE name and PAM']=[[val['BE name and editing window'].split(' editing window:')[0].replace('method:',''),
                             val['BE type and PAM'].split(' PAM:')[1]]]
    be_type=val['BE type and PAM'].split(' ')[0]                             
    cfg['BE type']=[[be_type.split('-')[0],be_type.split('-')[1]]]
    if not vals2 is None:
        cfg['PAM position']=[vals2['PAM position'].replace('stream','')]                    
        cfg['guide length']=[int(vals2['guide length'])]    

    window=val['BE name and editing window'].split(' editing window:')[1].replace('bp','')                         
    cfg['BE editing window']=[[int(window.split('-')[0]),int(window.split('-')[1])]]
    cfg['cores']=int(val['cores'])

    cfg['dinp']=val['mutation table']
    cfg['mutation_format']='nucleotide' if val['mutation_format nucleotide'] else 'aminoacid'
    cfg['reverse_mutations']=False if val['reverse_mutations create'] else True
    #advanced
#     cfg['dsubmap_preferred_path']=val['dsubmap_preferred_path']
#     cfg['mimetism_level']=val['mimetism_level']
#     keys=np.array([k for k in val if 'non-intermutable ' in str(k)])
#     buls=np.array([val[k] for k in keys])
#     non_intermutables=list(keys[buls])
#     non_intermutables=[s.replace('non-intermutable ','') for s in non_intermutables]
#     cfg['non_intermutables']=non_intermutables
    cfg['keep_mutation_nonsense']=True
    cfg['mutation_type']=None 
    cfg['make_control_pos']=val['positive controls']
    cfg['make_control_neg']=val['negative controls']
    cfg['mutations']='mutations'

#     deps=['samtools','bedtools','bwa',]
    cfg['gui']=True    
    if not vals2 is None:
        cfg['custom BE and PAM']=True
    else:
        cfg['custom BE and PAM']=False        
#     for dep in deps:
#         if val[dep]!='':
#             cfg[dep]=val[dep]
#         else:
#             cfg[dep]=dep            
    yaml.dump(cfg,open(val['cfgp'],'w'))
    return cfg
    # 

# gui io
def resetwinvals(win,vals,test=False): 
    for k in vals:
        try:
            # if not radio
            if not 'mutation_format ' in k:
                win.FindElement(k).Update(vals[k])                            
        except:
            if test:
                print(f'resetwinvals error with key={k}')
    radios= [['mutation_format aminoacid','mutation_format nucleotide'],
            ['reverse_mutations create','reverse_mutations remove']]               
    for options in radios:
        for opt, opt_ in zip(options,options[::-1]):
            if (opt in vals) and (opt_ in vals):
                if (not vals[opt]) and (not vals[opt_]):
                    win.FindElement(opt_).Update(False)
                    win.FindElement(opt).Update(False)
                    break 
                if vals[opt] and not vals[opt_]: 
                    win.FindElement(opt_).Update(vals[opt_])
                    win.FindElement(opt).Update(vals[opt])
                    break
    return win 

# def runcom(command, *args):      
#     try:      
#         sp = subprocess.Popen([command, *args], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)      
#         out, err = sp.communicate()      
#         if out:      
#             print(out.decode("utf-8"))      
#         if err:      
#             print(err.decode("utf-8"))      
#     except:      
#         pass
from beditor.lib.io_sys import runbashcmd

# gui aes 
def splitlist(l,n): return [l[i:i + n] for i in range(0, len(l), n)]
def h1(s,width=None,size=20,kws={}):return sg.Text(s, size=(int(len(s)*2) if width is None else width, 1), font=("Monospace", size),**kws)
def h2(s,width=15,size=15,kws={}):return sg.Text(f"  {s}", size=(int(len(s)*2) if width is None else width, 1), font=("Monospace", size),**kws)
def h3(s,width=None,size=10,kws={}):return sg.Text(f"    {s}", size=(int(len(s)*1.5) if width is None else width, 1), font=("Monospace", size),**kws)
def normal(s,width=None,size=10,kws={}):return sg.Text(f"{s}", size=(int(len(s)*1.5) if width is None else width, 1), font=("Monospace", size),**kws)
h2_font=("Monospace", 15)
h3_font=("Monospace", 10)
kws_frame={'font':h2_font,
'relief':'sunken'}
kws_button_big={'font':h2_font,}
kws_button_small={'font':h3_font,}

def get_layout(test=False):
    # gui aes
    sg.SetOptions(background_color='#FFFFFF',      
               text_element_background_color='#F8F8F8',      
               element_background_color='#F8F8F8',      
               scrollbar_color=None,      
               input_elements_background_color='#FFFFFF',      
               progress_meter_color = ('green', 'blue'),      
               button_color=('white','#6863FF'))

    win_width=80
    win = sg.Window('beditor', default_element_size=(win_width, 1),
    #                    resizable=True
                      )
    

    # beditorp=dirname(beditor.__file__)
    beditorp='.'

    dbepams=get_dbepams()

    dspecies=pd.read_table(f'{dirname(abspath(__file__))}/data/dspecies.tsv',keep_default_na=False)
    species=['Saccharomyces cerevisiae (R64-1-1)']+list(np.sort(dspecies['Scientific name (Ensembl Assembly)'].tolist()))

    layout_advanced_setting=[[sg.Button('load a configuration file', key='optional: load configuration file_',button_color=('black','white'))],    
                [sg.Frame('',
                [[sg.InputText(default_text='path to configuration file (.yaml)',key='cfginp'), 
                sg.FileBrowse(button_text='Browse',file_types=(('YAML file', '*.yml'),('YAML file', '*.yaml')),key='browse cfg',enable_events=False),
                sg.Button('Load', key='load cfg',disabled=False)],
                ],              
                tooltip='the yaml configuration file would fill the details of a beditor run.',
                key='optional: load configuration file',
                title_color='gray',
                visible=True if test else False,
                **kws_frame)],
#                 [sg.Button('specitf type of amino acid mutations', key='options for amino acid mutations_',button_color=('black','white'))],    
#                 [sg.Frame('',
#                           [[h3('mutation type'),
#                              sg.Checkbox('non-synonymous', key="mutation_type non-synonymous", default=True),
#                              sg.Checkbox('synonymous', key="mutation_type synonymous", default=True),
#                              sg.Checkbox('include non-sense', key="keep_mutation_nonsense", default=True)]],
#             #               tooltip='base editor/s to use.',
#                             visible=True if test else False,
#                             key='options for amino acid mutations',
#                           **kws_frame)],
#                 [sg.Button('infer mutated amino acid', key='optional: for infering amino acid mutations from wt_',
#                            button_color=('black','white'))],    
#                 [sg.Frame('',
#                         [[h3('substitution map'),
#                           sg.InputText(default_text='',key='dsubmap_preferred_path'),
#                           sg.FileBrowse(button_text='Browse',file_types=(('Tab separated values file (csv)', '*.csv'),)),
#                         sg.Text('', key='required substitution map',text_color='red',size=(25, 1)),
#                           ],
#                         [h3('mimetic mutations'),
#                         sg.InputCombo(['high','medium','low'],key='mimetism_level',default_value='mimetism_level',
#                         tooltip='mimetism based substitutions allow\nmutations to amino acids with similar properties as of wt.\nmimetism level (high: only the best one, [medium: best 5], low: best 10)',
#                         size=(15, 1)),
#                         sg.Text('', key='required mimetic mutations',text_color='red',size=(25, 1)),
#                         ],
#                         [h3('non-intermutable', width=20),]+[sg.Checkbox(aa, key=f"non-intermutable {aa}", default=False) for aa in aminoacids[:10]],
#                         [h3('amino acids     ', width=20),]+[sg.Checkbox(aa, key=f"non-intermutable {aa}", default=False) for aa in aminoacids[11:]]
#                         ],              
#                         tooltip='required (*) if amino acid mutation (to) is not provided.',
#                           title_color='gray',
#                         key='optional: for infering amino acid mutations from wt',
#                         visible=True if test else False,
#                           **kws_frame)],
                [sg.Button('design control gRNAs', key='optional: design control gRNAs_',button_color=('black','white'))],    
                [sg.Frame('',
                          [[sg.Checkbox('positive controls', key="positive controls", default=True),
                             sg.Checkbox('negative controls', key="negative controls", default=True),
                            ]],
            #               tooltip='base editor/s to use.',
                            visible=True if test else False,
                            key='optional: design control gRNAs',
                          **kws_frame)],                          
#                 [sg.Button('use custom dependencies', key='optional: dependencies paths_',button_color=('black','white'))],    
#                 [sg.Frame('',
#                     [[h3('bedtools',width=20),sg.InputText(default_text='bedtools',
#                                                            tooltip='bash command or path to the source',key='bedtools')],              
#                     [h3('samtools',width=20),sg.InputText(default_text='samtools',
#                                                           tooltip='bash command or path to the source',key='samtools')],              
#                     [h3('bwa',width=20),sg.InputText(default_text='bwa',
#                                                      tooltip='bash command or path to the source',key='bwa')],],
#                   tooltip='optional because the dependencies would be otherwise\ninstalled by beditor.',
#                           title_color='gray',
#                     visible=True if test else False,
#                             key='optional: dependencies paths',
#                           **kws_frame)]
                            ]
    # TAB01
    layout_configure=[
        # [h1('configure',width=40,kws={'tooltip':'this step would setup the parameters required for a beditor analysis.',})],
        [h2('genome',width=24,kws={'tooltip':'genome information of host organism.'}),
        sg.InputCombo(species,key='Species name (Ensembl assembly)',
                      default_value='Species name (Ensembl assembly)',
                tooltip='Species name (Ensembl assembly)',
                # enable_events=True,
                size=(60, 1)),     
#          sg.InputCombo([f"genome release={i}" for i in range(91,94,1)][::-1],tooltip='Ensembl genome release eg. 93',size=(17, 1),key='genomerelease')
        ],
        [sg.Frame('BE and PAM',
            [[h2('BE type and PAM',width=30,kws={'tooltip':'genome information of host organism.'}),
                sg.InputCombo(list(np.sort(dbepams['BE type and PAM'].unique())),
                        key='BE type and PAM',
                        disabled=False,
                        default_value='choose BE type and PAM',
                        enable_events=True,
                        size=(50, 1)),
                    ],
            [h2('BE name and editing window',width=30,kws={'tooltip':'genome information of host organism.'}),
                sg.InputCombo(list(np.sort(dbepams['BE name and editing window'].unique())),
                        key='BE name and editing window',
                        default_value='choose BE and editing window',
                        enable_events=True,
                        disabled=True,
                        size=(50, 1)),
                    ],
            [h2('',width=30),
            sg.Button('Custom',key='add_bepam',tooltip='Use custom new BE and PAM.'),
             sg.Text('', key='add_bepam print',text_color='green',size=(22, 1)),
             sg.Button('Clear', key='BE and PAM clear',disabled=True)],
            ],
            key='BE and PAM',
            tooltip='base editor/s to use.',
            **kws_frame)],    
        # [h1('step02 mutation info',width=40)],
        [h2('mutation table',width=24),
         sg.InputText(default_text='path to the tsv file',key='mutation table',size=(44, 1)), 
         sg.FileBrowse(button_text='Browse',file_types=(('Tab separated values file (tsv)', '*.tsv'),),
#             enable_events=True,key='browse din'
                      ),
        sg.Button('Load', key='load din',disabled=True),
        sg.Text('', key='error din',text_color='green',size=(25, 1)),
        ],    
        [h2('mutation mode',width=24),
         sg.Radio('create mutations', "mutation_do", key='reverse_mutations create'),
         sg.Radio('remove mutations', "mutation_do",key='reverse_mutations remove')],
        [h2('mutation format',width=24),
         sg.Radio('nucleotide', "mutation_format",key='mutation_format nucleotide',),
         sg.Radio('amino acid', "mutation_format",key='mutation_format aminoacid',)],
        [h2('assign # of cpus',width=24,kws={'tooltip':'number of cores/cpus to be used. higher is faster.'}),
         sg.Slider(range=(1, multiprocessing.cpu_count()-1 if multiprocessing.cpu_count()-1 !=0 else 1), orientation='h', size=(25, 5), default_value=int(multiprocessing.cpu_count()*0.5),key='cores'),
         sg.Checkbox('', key="calculate beditor scores", default=False),
         normal('beditor scores',size=15,width=30,kws={'tooltip':'celculate editability of gRNAs.'}),],
        [sg.Button('go to next step  '+' '*52, key='configuretorun',**kws_button_big),
        sg.Text('', key='configure error',text_color='red',size=(25, 1))],    

        ##OPTIONS
        [sg.Button('advanced settings', key='configure_advanced_',button_color=('black','white'),**kws_button_small),
                     sg.Button('Clear', key='clear all',disabled=False)],
        [sg.Frame('',
        layout_advanced_setting,
        # tooltip='optional because the dependencies would be otherwise\ninstalled by beditor.',
        title_color='gray',
        key='configure_advanced',
        visible=True if test else False, 
        **kws_frame)],
        ]  


    # TAB02
    layout_run = [
        # [h1('run beditor',width=40)],
    #     [sg.Text('_'  * width)],
        [normal('save configuration file',size=15)],
        [sg.InputText(default_text='path to save configuration file (.yml)',key='cfgp',size=(58,1)),
        sg.FileSaveAs(button_text='Browse',file_types=(('YAML file', '*.yml'),),size=(10,1),
#                       enable_events=True,
#                       key='browse cfgp'
                     ),
        sg.Button('Save', key='save cfgp',size=(10,1),disabled=True),
        sg.Text('', key='save cfgp error',text_color='red',size=(25, 1))
        ],
        [sg.Button('Run beditor'+' '*52, key='run beditor',**kws_button_big,disabled=True)],
        [sg.Image(f'{dirname(abspath(__file__))}/data/gui/guiload.gif',visible=True,key='guiload')],
        [sg.Text('', key='run beditor error',text_color='green',size=(50, 1))],    

        ]

    ## LAYOUT
    layout = [[sg.TabGroup([[sg.Tab('configure', layout_configure, tooltip='configure settings',key='configure',
                                    title_color='blue',font=("Monospace", 20),disabled=False),
                             sg.Tab('run beditor', layout_run, tooltip='run beditor',key='run',
                                    title_color='blue',font=("Monospace", 20),disabled=True if not test else False),
                            ]],
                            selected_title_color='blue')]]
    return layout

## POPUP
layout_addbepam = [
            [sg.Frame('BE',[[normal('name',width=25),sg.InputText(default_text='name',key='BE name')],
            [normal('nucleotide mutation',width=25,
            kws={'tooltip':'length of guide/sgrna sequence.'}),
            sg.Text('', key='error mutation',text_color='red',size=(25, 1))],
            [normal('from',width=25),
            sg.Radio('A', "wt nucleotide",key='fromA',default=False), 
            sg.Radio('T', "wt nucleotide",key='fromT',default=False),
            sg.Radio('G', "wt nucleotide",key='fromG',default=False),
            sg.Radio('C', "wt nucleotide",key='fromC',default=False)],
            [normal('to',width=25),
            sg.Radio('A', "mt nucleotide",key='toA',default=False), 
            sg.Radio('T', "mt nucleotide",key='toT',default=False),
            sg.Radio('G', "mt nucleotide",key='toG',default=False),
            sg.Radio('C', "mt nucleotide",key='toC',default=False)],
            [normal('editing window',width=25,
            kws={'tooltip':'length of guide/sgrna sequence.'}),
            sg.Text('', key='editing window error',text_color='red',size=(25, 1))],
            [normal('minimum distance from PAM',width=25,
            kws={'tooltip':'length of guide/sgrna sequence.'}),
            sg.Slider(range=(10, 24), orientation='h', size=(15, 5), default_value=13,key='editing window min'),],
            [normal('maximum distance from PAM',width=25,kws={'tooltip':'length of guide/sgrna sequence.'}),
            sg.Slider(range=(10, 24), orientation='h', size=(15, 5), default_value=20,key='editing window max'),],],**kws_frame)],
            [sg.Frame('PAM',[[normal('PAM',width=25),sg.InputText(default_text='PAM sequence',key='PAM'),
            sg.Text('', key='error',text_color='red',size=(25, 1))
            ],
            [normal('position',width=25),
            sg.InputCombo(['downstream','upstream'],default_value='upstream',tooltip='wrt position of gRNA',size=(15, 1),key='PAM position')],
            [normal('guide length',width=25,kws={'tooltip':'length of guide/sgrna sequence.'}),
            sg.Slider(range=(10, 24), orientation='h', size=(15, 5), default_value=20,key='guide length'),],
            ],**kws_frame)],  
            [sg.Button('add'),sg.Cancel()],  
            ]

    
def gui(test=False):
    layout=get_layout(test=test)
    win = sg.Window('beditor').Layout(layout)  
        
    win_addbepam_active=False  
    bulconfigure_advanced=False
    init=True
    while True:  
        # gottu be in while loop to capture first event
        ev1, vals1 = win.Read()
        if init:
            _ev1, _vals1 = ev1, vals1
            init=False
        if test:
            print(ev1)
            print(vals1)
        if vals1['mutation table']!='':
            win.FindElement('load din').Update(disabled=False)
        if vals1['cfgp']!='':
            win.FindElement('save cfgp').Update(disabled=False)
            
        if ev1 is None:  
            break  
        if ev1=='BE type and PAM':
            # print('event: BE type and PAM')
            dbepams=get_dbepams()
            dbepams_=dbepams.groupby('BE type and PAM').agg({'BE name and editing window':unique_dropna,}).reset_index()
            l=tuple(dbepams_.loc[dbepams_['BE type and PAM']==vals1['BE type and PAM'],'BE name and editing window'].sort_values().tolist()[0])
            win.FindElement('BE name and editing window').Update(disabled=False)
            win.FindElement('BE name and editing window').Update(values=l)
            win.FindElement('add_bepam').Update(disabled=True)                    
            win.FindElement('BE and PAM clear').Update(disabled=False)
            win.FindElement('BE name and editing window').Update(disabled=False)
            win=resetwinvals(win,vals1)
        elif ev1 == 'add_bepam'  and not win_addbepam_active:  
            win_add_bepam_active = True  
            win.Hide()  
            win_add_bepam = sg.Window('add new BE and PAM').Layout(layout_addbepam)  
            while True:            
                ev2, vals2 = win_add_bepam.Read() 
                if test:
                    print(ev2)
                if ev2 is None:  
                    win_add_bepam.Close()  
                    win_add_bepam_active = False  
                    win.UnHide()  
                    break
                elif ev2 == 'add':  
                    if any([vals2['fromA'] and vals2['toA'],vals2['fromT'] and vals2['toT'],
                            vals2['fromG'] and vals2['toG'],vals2['fromC'] and vals2['toC']
                           ]): 
                        if test:
                            print(ev2, vals2)
                        win_add_bepam=resetwinvals(win_add_bepam,vals2)
                        win_add_bepam.FindElement('error mutation').Update("* invalid")
                    elif not vals2['editing window min']<vals2['editing window max']:
                        win_add_bepam=resetwinvals(win_add_bepam,vals2)
                        win_add_bepam.FindElement('editing window error').Update("* invalid")
                    elif not isstrallowed(s=vals2['PAM'],form=f"^[{nts}]*$"):
                        win_add_bepam=resetwinvals(win_add_bepam,vals2)
                        win_add_bepam.FindElement('error').Update('* invalid nucleotide')
                    else:
                        # get the keys and print on the gui
                        if test:
                            print(vals2)
                        win_add_bepam.Close()
                        win_add_bepam_active = False
                        win.UnHide()
                        win.FindElement('add_bepam print').Update(f"{get_mutation(vals2)} {vals2['BE name']} {vals2['editing window min']}-{vals2['editing window max']}bp")
                        vals1['add_bepam print']=f"{get_mutation(vals2)} {vals2['BE name']} {vals2['editing window min']}-{vals2['editing window max']}bp"
                        win.FindElement('BE type and PAM').Update(values=[f"{get_mutation(vals2)} PAM:{vals2['PAM']}",])
                        win.FindElement('BE name and editing window').Update(values=[f"method:{vals2['BE name']} editing window:{int(vals2['editing window min'])}-{int(vals2['editing window max'])}bp",])
                        win.FindElement('BE type and PAM').Update(disabled=True)
                        win.FindElement('BE name and editing window').Update(disabled=True)
                        win.FindElement('BE and PAM clear').Update(disabled=False)
                        break
                elif ev2 == 'Cancel':
                    win_add_bepam.Close()  
                    win_add_bepam_active = False  
                    win.UnHide()  
                    break
            win=resetwinvals(win,vals1)
        elif ev1 == 'BE and PAM clear':
            win.FindElement('BE and PAM clear').Update(disabled=False)
            dbepams=get_dbepams()
            win.FindElement('BE type and PAM').Update(values=list(np.sort(dbepams['BE type and PAM'].unique())),disabled=False)
            win.FindElement('BE name and editing window').Update(disabled=False)
            win.FindElement('add_bepam').Update(disabled=False)
            win.FindElement('add_bepam print').Update("")
            win=resetwinvals(win,vals1)        
        elif ev1 == 'load cfg':
            if vals1['cfginp']!='' and vals1['cfginp']!='path to configuration file (.yaml)':
                cfg=yaml.load(open(vals1['cfginp'],'r'))
                listed_keys=['BEs','pams']
                del_keys=['max_subs_per_codon','mutations','chunksize']
                for k in listed_keys:
                    if not isinstance(cfg[k],list):
                        cfg[k]=[cfg[k]]
                for k in del_keys:
                    del cfg[k]
                win=loadcfginwinvals(win,cfg)
            win=resetwinvals(win,vals1)        
        elif ev1=='configure_advanced_':
            win=resetwinvals(win,vals1)
            win.FindElement('configure_advanced').Update(visible=False if bulconfigure_advanced else True)
            bulconfigure_advanced=False if bulconfigure_advanced else True
            win=resetwinvals(win,vals1)          
        elif ev1=='optional: load configuration file_':
            win.FindElement('optional: load configuration file').Update(visible=True)
            win=resetwinvals(win,vals1)          
        elif ev1=='options for amino acid mutations_':
            win.FindElement('options for amino acid mutations').Update(visible=True)
            win=resetwinvals(win,vals1)          
        elif ev1=='optional: for infering amino acid mutations from wt_':
            win.FindElement('optional: for infering amino acid mutations from wt').Update(visible=True) 
            win=resetwinvals(win,vals1)         
        elif ev1=='optional: dependencies paths_':
            win.FindElement('optional: dependencies paths').Update(visible=True) 
            win=resetwinvals(win,vals1)         
        elif ev1=='optional: or just save configuration file_':
            win.FindElement('optional: or just save configuration file').Update(visible=True)  
            win=resetwinvals(win,vals1)        
        elif ev1=='optional: design control gRNAs_':
            win.FindElement('optional: design control gRNAs').Update(visible=True)  
            win=resetwinvals(win,vals1)        
        elif ev1=='configuretorun':
            #check vals
            if vals1['mutation table']=='path to the tsv file':
                vals1['mutation table']=''    
            if win.FindElement('add_bepam print').DisplayText=='':
                keys=np.array(['mutation table','Species name (Ensembl assembly)',
                    'BE type and PAM','BE name and editing window'])
            else:
                keys=np.array(['mutation table','Species name (Ensembl assembly)'])
            buls=np.array([vals1[k]=='' for k in keys])            
            
            if len(keys[buls])==0:
                din=del_Unnamed(pd.read_table(vals1['mutation table'],keep_default_na=False))
                if (('genome coordinate' in din) and (vals1['mutation_format nucleotide'])) or (('transcript: id' in din) and (vals1['mutation_format aminoacid'])):                
                    win.FindElement('configure error').Update('run beditor',text_color='green')            
                    win.FindElement('run').Update(disabled=False)  
                else:             
                    win.FindElement('configure error').Update("invalid mutation format",text_color='red')
            else:             
                win.FindElement('configure error').Update(f"invalid {' and '.join(list(keys[buls]))}",text_color='red')
            win=resetwinvals(win,vals1)        
            ##TODO check columns and rename
        elif ev1=='load din':
    #         try:
            din=del_Unnamed(pd.read_table(vals1['mutation table'],keep_default_na=False))
            cols_din=['genome coordinate','nucleotide mutation','transcript: id','aminoacid: position','amino acid mutation']
            normalised2cols=dict(zip([normalisestr(c) for c in cols_din],cols_din))
            din=din.rename(columns={c:normalised2cols[normalisestr(c)] for c in din})
            if ('genome coordinate' in din) and (not 'transcript: id' in din):
                win.FindElement('mutation_format aminoacid').Update(value=False)
                win.FindElement('mutation_format nucleotide').Update(value=True)
                vals1['mutation_format aminoacid']=False
                vals1['mutation_format nucleotide']=True
            elif (not 'genome coordinate' in din) and ('transcript: id' in din):
                win.FindElement('mutation_format nucleotide').Update(value=False)
                win.FindElement('mutation_format aminoacid').Update(value=True)
                vals1['mutation_format aminoacid']=True
                vals1['mutation_format nucleotide']=False
            else:
                vals1['mutation_format aminoacid']=False
                vals1['mutation_format nucleotide']=False
            vals1['reverse_mutations create']=True
            vals1['reverse_mutations remove']=False   
            win.FindElement('error din').Update("loaded",text_color='green')
            win.FindElement('mutation table').Update(vals1['mutation table'])
            win=resetwinvals(win,vals1)
        elif ev1 == 'mutation_format nucleotide':
            vals1['mutation_format nucleotide']=True
            vals1['mutation_format aminoacid']=False
            win=resetwinvals(win,vals1)
        elif ev1 == 'mutation_format aminoacid':
            vals1['mutation_format nucleotide']=False
            vals1['mutation_format aminoacid']=True
            win=resetwinvals(win,vals1)
        elif ev1 == 'reverse_mutations create':
            vals1['reverse_mutations create']=True
            vals1['reverse_mutations remove']=False
            win=resetwinvals(win,vals1)
        elif ev1 == 'reverse_mutations remove':
            vals1['reverse_mutations create']=False
            vals1['reverse_mutations remove']=True
            win=resetwinvals(win,vals1)
        elif ev1=='clear all':
            win=resetwinvals(win,_vals1)

        elif ev1 == 'save cfgp':
            if vals1['cfgp']!='' or vals1['cfgp']!='path to save configuration file (.yml)':
                vals1['cfgp']= f"{vals1['cfgp']}.yml" if not vals1['cfgp'].endswith('.yml') else vals1['cfgp']
                if test:
                    yaml.dump(vals1,open(vals1['cfgp']+'_test.yml', 'w'), default_flow_style=False)
                if not 'vals2' in locals():
                    vals2=None
                cfg=guival2cfg(vals1,vals2)
                from beditor.pipeline import validcfg
                if validcfg(cfg):
                    yaml.dump(cfg,open(vals1['cfgp'], 'w'), default_flow_style=False)
                    win.FindElement('save cfgp error').Update(f"saved",text_color='green')
                else:
                    win.FindElement('save cfgp error').Update(f"error/s in configuration",text_color='red')
                win.FindElement('run beditor').Update(disabled=False)        
            win.FindElement('cfgp').Update(vals1['cfgp'])            
            win=resetwinvals(win,vals1)        
        elif ev1 == 'run beditor':      
            win.FindElement('guiload').UpdateAnimation(source=f'{dirname(abspath(__file__))}/data/gui/guiload.gif',
                                                       time_between_frames=0)
            win.FindElement('run beditor error').Update(f"running!",text_color='green')
            try:
                runbashcmd(f"source activate beditor; beditor --cfg {vals1['cfgp']}")
                win.FindElement('run beditor error').Update(f"finished processing!",text_color='green')
            except:
                win.FindElement('run beditor error').Update(f"errored! see command line",text_color='red')            
            #TODO create cfg and validate
            win=resetwinvals(win,vals1)        
        # elif ev1 == 'host':
        #     print(vals1['cfgoutp'])
        #     win.FindElement('mutation_format aminoacid').Update(True)
