import os

import htcondor


condordir = '/scratch365/atownse2/condor'
def submit_batch(executable, arguments, job_name):
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