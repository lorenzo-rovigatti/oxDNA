from sys import stderr

class LoggerSettings:
    def __init__(self, quiet:bool=0):
        self.quiet = quiet
    
    def set_quiet(self, value:bool):
        self.quiet = value

    def get_quiet(self):
        return self.quiet

logger_settings:LoggerSettings = LoggerSettings()

def log(message:str, end:str="\n", level:str|int="info"):
    """
        Log something to stdout.  Doesn't do anything fancy.
        
        Parameters:
            message(str) : The message to write
            level(str|int) : info(0)|warning(1) Whether this is informational or a warning.  For errors we just raise errors.
            quiet(bool) : if 1 only print warnings 
    """
    global logger_settings

    if isinstance(level, str):
        if level == "info":
            level = 0
        elif level == "warning":
            level = 1
        else:
            raise RuntimeError(f"{level} is not a recognized logging string")
    if logger_settings.get_quiet() and not level:
        return
    
    prefix:str = "WARNING" if level else "INFO"
    print(f"{prefix}: {message}", end=end, file=stderr)
    