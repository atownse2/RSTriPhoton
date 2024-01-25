import os

USER=os.environ['USER']

condordir = f'/scratch365/{USER}/RSTriPhoton/condor'

def check_condor_dir():
    if not os.path.isdir(condordir):
        os.makedirs(condordir)
        os.makedirs(f'{condordir}/out')
        os.makedirs(f'{condordir}/err')
        os.makedirs(f'{condordir}/log')

def submit_condor(executable, arguments, job_name):
    '''Submit a job to condor'''
    check_condor_dir()
    import htcondor
    make_events = htcondor.Submit({
        "executable": executable,
        "arguments": arguments,
        "output": f"{condordir}/out/{job_name}.out",
        "error" : f"{condordir}/err/{job_name}.err",
        "log"   : f"{condordir}/log/{job_name}.log",              
        "request_cpus": "1",
        "request_memory": "128MB",
        "request_disk": "128MB",
    })
    print(f"Submitting job with name : {job_name}")
    schedd = htcondor.Schedd()
    submit_result = schedd.submit(make_events)
    return submit_result

def submit_interactive(executable, arguments):
    command = f'{executable} {arguments}'
    os.system(command)
    return 0