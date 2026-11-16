"""Execute the standalone quantitative notebook in the Pixi kernel.

Uses installed jupyter_client/ipykernel, writes execution records under outputs,
and leaves the source notebook untouched.
"""
import argparse
import json
from pathlib import Path
import queue
import sys
import time

from jupyter_client import KernelManager


def run(output,notebook_path='notebooks/quantitative_headless_simulation.ipynb',tag=None):
    root=Path(__file__).resolve().parents[1]
    notebook=json.loads((root/notebook_path).read_text())
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    km=KernelManager(kernel_name='python3')
    km.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
    km.start_kernel(cwd=str(root))
    client=km.client();client.start_channels();records=[]
    try:
        client.wait_for_ready(timeout=60)
        for index,cell in enumerate(notebook['cells']):
            if cell['cell_type']!='code' or (tag and tag not in cell.get('metadata',{}).get('tags',[])):continue
            msgid=client.execute(''.join(cell['source']))
            deadline=time.monotonic()+120;outputs=[]
            while time.monotonic()<deadline:
                try:message=client.get_iopub_msg(timeout=1)
                except queue.Empty:continue
                if message.get('parent_header',{}).get('msg_id')!=msgid:continue
                kind=message['msg_type'];content=message['content']
                if kind=='error':raise RuntimeError('\n'.join(content['traceback']))
                if kind=='stream':outputs.append(content['text'])
                if kind=='status' and content['execution_state']=='idle':break
            else:raise TimeoutError(f'notebook cell {index}')
            reply=client.get_shell_msg(timeout=10)
            if reply['content']['status']!='ok':raise RuntimeError(reply['content'])
            records.append({'cell':index,'status':'passed','text_output':outputs})
            print(f'Cell {index}: passed',flush=True)
    finally:
        client.stop_channels();km.shutdown_kernel(now=True)
    if not records:raise ValueError('no matching notebook code cells')
    result={'scope':f'{notebook_path}; tag={tag!r}; only recorded cells executed',
            'python':sys.executable,'cells':records,'fitting':'not invoked'}
    (output/'notebook_execution.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--notebook',default='notebooks/quantitative_headless_simulation.ipynb')
    parser.add_argument('--tag',help='Execute code cells with this notebook tag only')
    args=parser.parse_args()
    run(args.output,args.notebook,args.tag)
