import logging
import sys
import os

class mapDamageLogger(logging.Logger):
    """
    Logger specific for mapDamage. 
    """
    def specialize(self, fo):
        """
        Specialize the Logger 
        """
        file_name=os.path.join(fo,"Runtime_log.txt")
        #Create a file handle
        file_handler = logging.FileHandler(file_name)
        file_handler.setLevel(0) #Set the level of the output
        #Formatter of the string 
        format_handler = logging.Formatter('%(asctime)s\t\t%(levelname)s\t\t%(message)s')
        file_handler.setFormatter(format_handler)
        self.addHandler(file_handler)
    def commandline(self):
        """
        Write the command line in the log file
        """
        self.info("Started with the command: "+" ".join(sys.argv))
    def finish(self):
        """
        Completed the computation
        """
        self.info("Successful run")

