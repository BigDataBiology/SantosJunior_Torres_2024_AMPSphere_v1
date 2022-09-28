import sys, select

def timeout_input(prompt, timeout=3, default=""):
    print(prompt, end=': ', flush=True)
    inputs, outputs, errors = select.select([sys.stdin], [], [], timeout)
    print()
    return (0, sys.stdin.readline().strip()) if inputs else (-1, default)

