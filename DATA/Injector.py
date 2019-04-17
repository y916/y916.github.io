#!/usr/bin/python
# -*- coding: utf-8 -*-

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import os

class Injector(object):
    def __init__(self,Density0,Material0,WP0,MDS,):
