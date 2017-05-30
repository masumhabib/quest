import os, sys
import configparser
import constants

WORK_DIR        = "work"
REPORT_DIR      = "report"

class TestcaseConfig(object):
    def __init__(self, conf_file):
        self.path = os.path.dirname(os.path.realpath(conf_file))
        self.filename = os.path.basename(conf_file)
        self.cfg = configparser.ConfigParser()
        self.cfg.read(conf_file)

    def save(self):
        self.cfg.write(self.path + "/" + self.filename)

    def input_file():
        return self.cfg['input_file']


