"""EGIS Python Utilities

Utility functions for use with arcpy to ease arcpy script tool development.
This module is designed to be a general purpose helper module
for ArcGIS python scripting.

This module supports both ArcGIS 9.x (arcgisscripting) and 10.x (arcpy).
It been implemented to support Python 2.x and 3.x.

Example scripts that use this module are included as part of the module

>>> print egis.TEMPLATE   (arcpy)
>>> print egis.TEMPLATE9  (arcgisscripting)

Author: Curtis Price, curtvprice@gmail.com
"""
# python 3.x compatibility check
from __future__ import print_function, unicode_literals, absolute_import

import sys
import os
import time
import traceback
import arcgisscripting
import arcpy

# version
Version = r"EGIS Python Utilities 11/14/2014 4:39:26 PM"

# script folder
__file__ = os.path.realpath(__file__)
ScriptPath = os.path.dirname(__file__)

# ArcGIS 10 distribution EGIS Toolbox path
# here:    ArcToolbox\Toolboxes\USGS_EGISTools\scripts
# toolbox: ArcToolbox\Toolboxes\USGS EGIS Tools.tbx
Toolbox = os.path.normpath(
    os.path.join(ScriptPath,"../..","USGS EGIS Tools.tbx"))
if not os.path.exists(Toolbox):
    # ArcGIS 9.x distribution EGIS toolbox path
    # here:    ArcToolBox\Scripts_egis
    # toolbox: ArcToolBox\Toolboxes\USGSEGISTools.tbx
    Toolbox = os.path.normpath(
        os.path.join(ScriptPath,"../Toolboxes","USGSEGISTools.tbx"))
    if not os.path.exists(Toolbox):
        Toolbox = None

# initialize GP variable
gp = None

# Set up timer function for GPMsg("t")
# Thanks to: Matthew Collier, USGS student, mcollier@usgs.gov
import time
global dt
t0 = time.clock()      # for execution time reporting at end of script
dt = [0,t0]            # set up a time object to report on differences
def timeDiff(T=None):  # define function to juggle and print time differences
    """Timer function for use in GPMsg()"""
    if not T: T=dt
    T[0] = T[1]
    T[1] = time.clock()
    return T[1] - T[0]


def getGP(GPVersion="9.3",strExt=None):
    """Instantiates a geoprocessor object.

    Arguments

      GPVersion - GP version, specified as number or string
        valid values:
        9.2 - http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?id=729&pid=727&topicname=Creating_the_geoprocessor_object
        9.3 (default) - http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?id=916&pid=914&topicname=Creating_the_geoprocessor_object
        10.0, 10.1 (not recommended -- use arcpy instead)

      strExt - ";" - delimited list of extensions

      NOTE: ESRI recommends caution mixing gp and arcpy in  multi-threads:
        FAQ:  Can multiple geoprocessor objects be created and used within an application?
        http://resources.arcgis.com/content/kbase?fa=articleShow&d=38286

    Examples

      gp = getGP()              # create a 9.3 gp
      gp = getGP(9.2)          # create a  9.2 gp
      gp = getGP(9.3,"spatial") # create a 9.3 GP and check out spatial analyst

      Note: This function is not needed for arcpy scripts.
    """
    global gp
    try:
        GPVersion = float(GPVersion)
        if int((10 * GPVersion) + 0.1) not in (92,93,100,101):
            GPMsg("w","getGP: Geoprocessor version %.1f" % GPVersion + \
                  " is not supported. Using 9.3 geoprocessor.")
            raise
    except:
        GPVersion = 9.3
    gp = arcgisscripting.create(GPVersion)

    if strExt:
        for ex in strExt.split(";"):
            try:
                strStatus = CheckOutExtension(ex)
                if strStatus != "CheckedOut":
                    GPMsg("w",ex)
            except Exception as xmsg:
                GPMsg("w",str(xmsg))
                pass
    return gp

def SetProduct(strProduct="ArcInfo"):
    """Run gp.SetProduct,  returning an error is not successful

    example

      SetProduct("ArcInfo")
    """
    global gp
    if not gp: gp = getGP()

    Result = None
    try:
        Result = gp.SetProduct(strProduct)
        if Result not in ["CheckedOut","AlreadyInitialized"]: raise Exception
    except:
        strMsg = "%s license not available." % strProduct
        raise MsgError(strMsg)
    return Result

def ArcGISProfile(profilePath=""):
    """Returns the ArcGIS Desktop user profile pathname

    The ArcGIS User profile includes various files such as
    user styles, toolboxes, and templates.

    example

      >>> print(ArcGISProfile("ArcToolbox"))
      C:\Documents and Settings\cprice\Application Data\ESRI\Desktop10.0\ArcToolbox\

    """
    # modified to support 10.1 4/17/2012

    global gp
    if not gp: gp = getGP()

    # get application data from Windows env variable
    APPDATA = os.environ["APPDATA"]
    # Check which version we've got.
    ArcVersion =  gp.GetInstallInfo("Desktop")['Version']
    if ArcVersion.find("10") == -1:
        # user settings path for ArcGIS 9.x
        path = os.path.join(APPDATA,"ESRI")
    else:
        # user settings path for ArcGIS 10.x
        path = os.path.join(APPDATA,"ESRI/Desktop%s/" % ArcVersion)
    path = os.path.realpath(os.path.join(path,profilePath))

    if not gp.exists(path): GPMsg("w","%s does not exist" % path)

    return path


from decimal import *
from math import log10, floor
def USGSRound(x, rval=0, rtype="DECIMALS"):
    """Return 'USGS-Rounded' numeric value (as Decimal)

    arguments

       x      string or numeric value
       rval   rounding value (0,1,2,3, ... 20)
       rtype  DECIMALS - number of decimals (default)
              SIGFIGS - number of significant figures (rval > 0)

    example:

        >>> for k in [-.5, 0, 99.45, 99.55, 253.8]:
        ...           print(k, " rounds to one place as ", USGSRound(k,1))
        ...           print(k, " rounds to one sigfig as ", USGSRound(k,1,"SIG"))
        ...
        -0.5  rounds to one place as  -0.5
        -0.5  rounds to one sigfig as  -0.5
        0  rounds to one place as  0
        0  rounds to one sigfig as  0
        99.45  rounds to one place as  99.4
        99.45  rounds to one sigfig as  100
        99.55  rounds to one place as  99.6
        99.55  rounds to one sigfig as  100
        253.8  rounds to one place as  253.8
        253.8  rounds to one sigfig as  300
    """
    if float(x) == 0.0:
        return Decimal('0')

    xstr = str(x)
    try:
        rval = int(rval)
        if rval > 20: raise
    except:
        raise Exception("invalid rval: %s" % rval)
    if rtype.lower()[:1] == "d":
        # number of decimals
        # convert number of decimals to quantize value (e.g. 2 -> .01)
        xrval = str(10 ** (-rval))
        xd = Decimal(xstr).quantize(Decimal(xrval), ROUND_HALF_EVEN)
    else:
        # number of sig figs
        if rval <= 0:
            raise Exception("invalid rval: %s" % rval)
        xflt = abs(float(x))
        mag = int(floor(log10(xflt)))
        fact = 10 ** (mag - rval + 1)
        xd = float(USGSRound(xflt / fact)) * fact
        xd = USGSRound(xd, max(0, rval - mag - 1))
        if float(x) < 0: xd = xd * Decimal('-1')
    return xd


def ScratchName(strPrefix="",strSuffix="",strDataType="",strWorkspace=""):
    """Return a scratch name with a prefix and suffix.

    This method is a wrapper of the gp.CreateScratchName
    method with a few enhancements.

    arguments

      prefix - prefix for scratch name
      suffix - suffix added to scratch name (default="")
      dataType - includes CreateScratchName options:

    notes

    1) If a workspace is not provided, the scratch workspace
    is determined by picking the first valid workspace of:
    scratch workspace, current workspace, or Windows TEMP.
    If the currently set scratch or current workspace is a
    personal or file geodatabase, and the request scratch name
    is for a file-based path (folder, grid, shapefile, etc),
    the scratchname will be generated for the folder above
    the geodatabase. ("path/xxx.gdb -> path/xxshp.shp")

    2) The scratch name returned is verified against the workspace
    to ensure that the pathname does not have a conflict, for example
    if creating a folder, it checks to ensure there is not an existing
    file with the same name. For coverages and grids, a check is made
    of the INFO directory for name conflicts.

    3) If the suffix includes an file extension, and the datatype
    supports extensions (for example, folders and coverages do not),
    the extension is returned.
      >>> gp.CreateScratchName("",".img","raster","e:/work")
      u'e:/work\\xx0'
      >>> nact.ScratchName("",".img","raster","e:/work")
      u'e:/work\\xx0.img'

    """
    global gp
    if not gp: gp = getGP()

    strDT = strDataType.lower()

    # Check out prefixes and suffixes
    if strPrefix == "": strPrefix = "x"
    pp = os.path.splitext(strPrefix)
    if pp[1] and strSuffix == "":
        # split prefix into prefix + suffix
        strPrefix = pp[0]
        strSuffix = pp[1]
    strPre = strPrefix.lower()
    strSuf = strSuffix.lower()
    # don't let it start with a number
    try:
        if strPre[0].isdigit(): strPre = "x" + strPre
    except:
        pass

    # make sure our scratch workspace is a folder if required
    strDT = strDataType.lower()

    # add .dbf suffix for dbase files
    if strSuf == ".dbf": strDT = "dbase"
    if strDT == "dbase" and strSuf.find(".dbf") == -1:
        strSuf += ".dbf"

    if strDT in ["coverage","folder","shapefile","arcinfotable",
                 "workspace","dbase","grid","tin"]:
        strScratchWS = ScratchWS(strWorkspace,"Folder")
        # change strDataType for gp.CreateScratchName syntax
        if strDT == "grid": strDataType = "coverage"
        elif strDT == "tin": strDataType = "folder"
    else:
        strScratchWS = ScratchWS(strWorkspace)

    # validate prefix, suffix names against workspace
    vff = strPre + strSuf
    try:
        if vff[0] == "." or vff[0].isdigit():
            strPre = "x"
            strSuf = vff[0]
    except:
        pass

    # validate prefix
    strPre = gp.ValidateTableName(strPre,strScratchWS)

    # loop until we're SURE we have an available pathname
    InfoList = None # Initialize - in case we need it (see below)
    GotValidScratchName = False    # we'll loop until this is true
    TryNum = 0 # incrementor
    while not GotValidScratchName:
        # set incremented scratch name prefix
        if TryNum == 0:
            strPre1 = strPre
        else:
            strPre1 = strPre + str(TryNum)
        # invoke gp.CreateScratchName method
        strScratchName = \
                       gp.CreateScratchName(strPre1,strSuf,strDataType,strScratchWS)
        # add back extension if gp.CreateScratchName dropped it
        # and "." is allowed in name
        SplitPath = os.path.splitext(strScratchName)
        if  SplitPath[1] != strSuf and strDT not in \
            ["coverage","grid","tin","folder"]:
            strScratchName += strSuf
            SplitPath = os.path.splitext(strScratchName)
        # Is this a valid scratch name?
        GotValidScratchName = True  # innocent until proven guilty!
        # does it exist already?
        if gp.Exists(strScratchName):
            GotValidScratchName = False
        # if shapefile, does .dbf
        elif strDataType == "shapefile" or SplitPath[1] == ".shp":
            # check for .dbf named same as proposed .shp file
            if gp.Exists(SplitPath[0] + ".dbf"):
                GotValidScratchName = False
        elif strDT in ["coverage","grid","tin","folder"]:
            # validating the name
            strScratchNameRoot = os.path.basename(strScratchName)
            if strDataType == "folder":
                strScratchNameRoot = strScratchNameRoot[:12]
            gp.ValidateTableName(strScratchNameRoot,os.path.dirname(strScratchName))
            strScratchName = os.path.join(\
                os.path.dirname(strScratchName),strScratchNameRoot)
            if strDT != "folder":
                # check for name conflicts in INFO tables "<DSName>.xxx"
                if not InfoList:
                    envWS = gp.Workspace # remember workspace
                    gp.Workspace = strScratchWS
                    InfoList = gp.ListTables("*","ArcInfoTable")
                    InfoList = [os.path.splitext(k)[0].lower() for k in InfoList]
                    gp.Workspace = envWS # restore
                try:
                    # if we get an error, there was a match!
                    idx = InfoList.index(strScratchNameRoot.lower())
                    GotValidScratchName = False
                except:
                    pass

        TryNum += 1  # increment, for next time around
    return os.path.normpath(strScratchName)

def ScratchFolder():
    """env.scratchFolder - for all versions of gp/arcpy
    1) If supported, return arcpy.env.scratchFolder (ArcGIS 10.1 or later)
    2) Return a value based on ScratchWorkspace (ArcGIS 9.3 - 10.0)
        a) Folder: scratchWorkspace
        b) GDB: scratchWorkspace/../scratch
        c) Not set or invalid: TEMP/scratch
    """
    global gp
    if not gp: gp = getGP()

    try:
        sw = gp.scratchWorkspace
        if sw:
            # check for invalid scratchWorkspace path
            if not gp.Exists(sw):
                GPMsg("w", "ID 873 {0}|\"{1}\"".format(
                    "Scratch Workspace invalid", sw))
                raise
        sw = gp.scratchFolder
    except:
        try:
            sw = gp.scratchWorkspace
            try:
                swType = gp.Describe(sw).dataType
                if swType == "Folder":
                    pass
                elif swType == "Workspace":
                    pth = os.path.dirname(sw)
                    sw = os.path.join(pth,"scratch")
                else:
                    raise
            except:
                sw = os.path.join(os.environ["TEMP"],"scratch")
        finally:
            sw = os.path.realpath(sw)
            if not gp.Exists(sw): os.mkdir(sw)
    return sw

def ScratchGDB():
    """env.scratchGDB for all versions of arcpy/gp

    1) If supported, return arcpy.env.scratchGDB (ArcGIS 10.1 or later)
    2) Return a value based on ScratchWorkspace (ArcGIS 9.3 - 10.0)
        a) GDB: scratchWorkspace
        b) Folder: scratchWorkspace/scratch.gdb
        c) Not set or invalid: TEMP/scratch.gdb
    """
    global gp
    if not gp: gp = getGP()

    try:
        sw = gp.scratchWorkspace
        if sw:
            # check for invalid scratchWorkspace path
            if not gp.Exists(sw):
                GPMsg("w", "ID 873 {0}|\"{1}\"".format(
                    "Scratch Workspace invalid", sw))
                raise
        sw = gp.scratchFolder
    except:
        try:
            sw = gp.scratchWorkspace
            swType = gp.Describe(sw).dataType
            if swType == "Workspace":
                pass
            elif swType == "Folder":
                sw = os.path.join(sw,"scratch.gdb")
        except:
            sw = os.path.join(os.environ["TEMP"],"scratch.gdb")
        finally:
            sw = os.path.realpath(sw)
            if not gp.Exists(sw):
                gp.CreateFileGDB_management(os.path.dirname(sw),
                                               os.path.basename(sw))
    return sw


def BuildQuery(table, field, operator="IS NOT NULL", value="", quote=None):
    """Simple SQL query expression generator

    arguments

      table     input feature class or table view
      field     field name (string)
                (The table and field names are not validated
      operator  SQL operator ("=","<>", etc)
                "IS NULL" and "IS NOT NULL" are supported
                if IN and IS NULL are used, value will not be quoted
      value     query value (optional)
      quote     Boolean - force quote or unquote value.
                if None (default) determine quoting from value

    examples

    BuildQuery("eggs.dbf", "avail", "=", "spam")
    '"AVAIL" = \'spam\''
    BuildQuery("spam.mdb/eggs", "ANSWER", ">", 42)
    '[ANSWER] > 42'
    BuildQuery("eggs.mdb", "ANSWER", "IN", "('SPAM', 'SPAM')")
    "[ANSWER] IN ('SPAM', 'SPAM')"
    BuildQuery("spam.dbf", "EGGS", "=", "{0}")
    '"EGGS" = \'{0}\''
    BuildQuery("spam.dat", "EGGS", "=", "{0}", False)
    '"EGGS" = {0}'
    """
    global gp
    if not gp: gp = getGP()

    try:

        # Use correct field delimeter (depends on table type)
        # and capitalize field name (ARCGIS SQL style)
        field = gp.AddFieldDelimiters(table, field.upper())

        # uppercase operator
        operator = operator.upper()

        # tweak string for selection value
        if quote == None:
            # determine quote from content
            if (isinstance(value, basestring)) and \
               operator not in ("IN","NOT IN"):
                quote = True
        if operator in ("IS NULL","IS NOT NULL"):
            # drop value when not specified (used for IS NULL, etc.)
            value = ""
            quote = False

        # build expression
        if quote:
            # add quotes around string values
            value = "'{0}'".format(value)

        sql = "{0} {1} {2}".format(field, operator, value)
        return sql.strip()

    except:
        raise MsgError("Invalid input")




# Functions for decimal degree / degrees, minutes seconds conversion
# Author: Curtis Price, http://profile.usgs.gov/cprice
# Disclaimer: Not approved by USGS -- Provisional, subject to revision.

def DDToDMS(dd, returnType="TEXT", numDec=3, Flip="NOFLIP"):
    """Convert a decimal degree coordinate pair to degrees, minutes seconds.

    Arguments

      dd - decimal-degree value (number)

      ReturnString -
         "TEXT": return a degrees-minutes-seconds string (Default)
         "NUM" return a 3-tuple of degrees,minutes,seconds numbers
         "TEXTNUM": return a degrees-minutes-seconds string (no spaces)

      NumDec -If returning strings, use this precision for the
              seconds (default: 3)

      If provided as negative a fixed format will be used
      (trailing zeros will be included).


      Flip -
            "NOFLIP" : Use the DD value as-is. (default)
            "FLIP"   :  Flip the sign by mulitplying the decimal degree
                   value by -1 before converting.
                   (Often used for longitudes in the western hemisphere)

    Examples

      >>> DDToDMS(-74.234)
      '-74 14 2.4'
      >>> DDToDMS(-74.234,"TEXT",-3)
      '-74 14 02.400'
      >>> DDToDMS(-74.234,"TEXT",3)
      '-74 14 02.4'
      >>> DDToDMS(-74.234,"TEXT",-3)
      '-74 14 02.400'
      >>> DDToDMS(-74.234,"TEXTNUM",-3)
      '-741402.400'
      >>> DDToDMS(-74.234,"NUM")
      (-74, 14, 2.3999999999999999)

      Author: Curtis Price, http://profile.usgs.gov/cprice
      (Not approved by USGS -- Provisional, subject to revision.)
    """
    try:
        returnType = str(returnType).upper()
        if returnType not in ["TEXT","NUM","TEXTNUM"]:
            raise Exception ("Option invalid: {0}".format(returnType))
        # Calculate deg,min,sec
        dd1 = abs(float(dd))
        xdeg = int(dd1)
        xmin = int(((dd1 - xdeg ) * 60.0) + 0.0001)
        # round seconds
        xsec = round(abs(( dd1 - xdeg - xmin / 60.0 ) * 3600),abs(numDec))
        if float(dd) < 0: xdeg = xdeg * -1
        # flip sign if requested
        if str(Flip).upper() in ["FLIP","TRUE"]:
            xdeg = xdeg * -1

        # return DMS value
        if returnType[:4] == "TEXT":
            # format to a string representation'
            # whole DMS first
            DMS = "%i %02i %02i" % (xdeg,xmin,int(xsec))
            if numDec != 0:
                if numDec > 0:
                    fcode = "g"
                else:
                    fcode = "f"
                fmt = " %%.%s%s" % (abs(numDec),fcode)
                # add .nnn
                DMS += (fmt % (xsec % 1))[2:]
            if returnType == "TEXTNUM":
                DMS = DMS.replace(" ","")
            return DMS
        else:
            # return numeric representation
            return xdeg,xmin,xsec
    except:

        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        print(tbinfo)
        raise Exception("{0}\nInvalid input".format(tbinfo))


def DMSToDD(dmsValue,Flip="NOFLIP"):
    """Convert a degrees,minutes,second string to decimal degrees

    Arguments

      dmsValue  Degrees-minutes-seconds value

                Valid formats include:

                "dd mm ss.ss" (text)
                ddmmss.ss (number or text)

      flip      Flip coordinates:

            "NOFLIP" : Use the DD value as-is. (default)
            "FLIP"   :  Flip the sign by mulitplying the decimal degree
                   value by -1 before converting.
                   (Often used for longitudes in the western hemisphere)


    Examples

        >>> DMSToDD(743000,"flip")
        -74.5
        >>> DMSToDD("74 30 00")
        74.5
        >>> DMSToDD(742000,"flip")
        -74.33333333333332

    """
    # Author: Curtis Price, http://profile.usgs.gov/cprice
    # Disclaimer: Not approved by USGS. (Provisional, subject to revision.)
    try:
        # convert dd mm ss.ss to dddmmss.sss
        Splits = [float(k) for k in str(dmsValue).split()]
        if len(Splits) == 3:
            xdms = "%i%02i%02i" % (Splits[0],Splits[1],int(Splits[2])) +\
                 ("%g" % (Splits[2] % 1))[2:]
        elif len(Splits) == 1:
        # if we didn't have 3-tuples, just remove any spaces
            xdms = str(dmsValue).replace(" ","")
        else:
            raise

        # okay, we have a string in the format: "ddmmss.sss"

        # calculate degrees, minutes, seconds
        xdeg = int(float(xdms) / 1e4)
        if xdeg < 0: sign = -1
        else: sign = 1
        # find degree value in xdms string
        # with that, we can get the minutes and seconds by position
        minpos = xdms.find(str(xdeg)) + len(str(xdeg))
        xmin = float(xdms[minpos:minpos+2])
        xsec = float(xdms[minpos+2:])

        # DMS triplet to decimal degrees
        dd = abs(xdeg) + xmin / 60.0 + xsec / 3600.0
        dd = dd * sign

        # flip sign if requested
        if str(Flip).lower() in ["flip","true"]:
            dd = dd * -1

        return dd

    except Exception as xmsg:
        raise Exception("Invalid input")

def HasIndex(tableView,fieldName):
    """Is a field indexed? (returns boolean)"""
    global gp
    if not gp: gp = getGP()

    tblPath = gp.Describe(tableView).CatalogPath
    fieldName = fieldName.upper()
    for i in gp.Describe(tblPath).Indexes:
        for field in i.Fields:
            if field.name.upper() == fieldName:
                return True
    return False

def AddIndex(tview,field,name=None,uniq=None,ascend=None):
    """Robustly try to add attribute index regardless of table type"""
    global gp
    if not gp: gp = getGP()

    try:
        gp.AddIndex_management(tview,field,name,uniq,ascend)
    except:
        try:
            # probably a shapefile, last 3 arguments not supported
            gp.AddIndex_management(tview,field)
        except:
            pass


def ScratchWS(strWS="",strWSType=""):
    """Validate a scratch workspace path

    arguments

    strWS (optional) - existing workspace path

    strWSType (optional)
      - "Folder" - always return a folder path
      - "Geodatabase" - always return a geodatabase
      - "in_memory" - always return the "in_memory" workspace

    notes

    The first valid path in this list is returned:
      1) strWS (if strWS is a valid path)
      2) gp.ScratchWorkspace
      3) gp.Workspace
      4) If a gdb is found (.gdb/.mdb )and strWSType == "Folder":
           folder above workpace ("xxx.gdb/..")
      5) Windows TEMP folder

    """
    global gp
    if not gp: gp = getGP()

    # parse workspace type
    strWST = strWSType.lower()
    if strWST.find("f") != -1: strWSType = "Folder"
    elif strWST.find("g") != -1: strWSType = "Geodatabase"

    # in_memory is a special case - just return it
    elif strWST.find("m") != -1 or strWS.lower() == "in_memory":
        return "in_memory"

    # find a valid scratch workspace path
    if not gp.Exists(strWS):
        strWS = gp.ScratchWorkspace
    if not gp.Exists(strWS):
        strWS = gp.Workspace
    if not gp.Exists(strWS):
        strWS = gp.GetSystemEnvironment("TEMP")
    # Check the path we have found to make sure it aggress
    # with the strWSType argument.
    WSExt = os.path.splitext(strWS)[1]
    WSExt = WSExt.replace(".","").lower() # ".gdb" -> "gdb"
    if strWSType == "Folder":
        if WSExt in ["mdb","gdb"]:
            strWS = os.path.dirname(strWS)
        if not gp.Exists(strWS):
            strWS = os.environ["TEMP"]
    elif strWSType == "Geodatabase":
        if WSExt not in ["mdb","gdb"]:
            raise Exception("%s is not a geodatabase workspace" % strWS)
    return strWS


def GetExtent(Dataset=None,Grace=0.0):
    """Returns a geoprocessing extent (as a string)

    arguments

    Dataset (Feature layer or geodataset) -
      * If not specified, the current GP extent is used.
      * If the feature layer has an active selection,
        only those features are used to determine the extent.

    Grace (number) -
      length to expand the extent (in extent's units)

    notes

    This function will not alter the current
    extent, but you can use it do so:

    gp.Extent = egis.GetExtent("e:/work/mygrid")
    gp.Extent = egis.GetExtent("polygon layer")

    Arc 10.0
    arcpy.env.extent = egis.GetExtent("e:/work/mygrid")

    """
    global gp
    if not gp: gp = getGP()

    try:
        lstExt = None  # initialize variable
        d = gp.Describe(Dataset)  # if dataset not valid, go to except
        Ex = gp.Describe(d.CatalogPath).Extent # returns an extent object
        lstExt = [Ex.XMin, Ex.YMin, Ex.XMax,Ex.YMax]
        if d.DataType != "FeatureLayer": raise Exception # our work is done
        # if this IS a layer - make sure there is a selected set
        numSel = int(gp.GetCount(Dataset).GetOutput(0))
        numRow = int(gp.GetCount(d.CatalogPath).GetOutput(0))
        if numSel == numRow: raise Exception # nothing selected, our work is done
        # Okay, this is a feature layer with selected features
        # Find the extent of selected features
        Rows = gp.SearchCursor(lyrFeatures)
        Row = Rows.Next()
        Ex = Row.Shape.Extent
        lstExt = [Ex.XMin, Ex.YMin, Ex.XMax,Ex.YMax]
        Row = Rows.Next()
        while Row:
            Ex = Row.Extent
            if lstExt[0] > Ex.XMin: lstExt[0] = Ex.XMin
            if lstExt[1] > Ex.YMin: lstExt[1] = Ex.YMin
            if lstExt[2] < Ex.XMax: lstExt[2] = Ex.XMax
            if lstExt[3] < Ex.YMax: lstExt[3] = Ex.YMax
            Row = Rows.Next()
            del Row, Rows
    except Exception as xmsg:
        if lstExt == None:
            # this is not a data set - use current extent
            try:
                strExtent = gp.Extent
                lstExt = [float(k) for k in strExtent.split()]
            except:
                # if we couldn't parse it, the extent is something like "MINOF"
                if Grace != 0:
                    # can't add a grace area to something like "MINOF", sorry
                    raise Exception(
                          "Can't add grace area to extent \"%s\"" % strExtent)
                return strExtent
        if not gp.Exists(d.CatalogPath):
            # something unexpected went wrong
            raise Exception(xmsg)

    # Add "grace area" if requested
    try:
        Grace = float(Grace)
        if Grace != 0:
            lstExt[0] -= Grace
            lstExt[1] -= Grace
            lstExt[2] += Grace
            lstExt[3] += Grace
    except:
        raise Exception("Invalid value for Grace")
    lstExt = [str(ee) for ee in lstExt]
    strExtent = "%s %s %s %s" % tuple(lstExt)

    # round numbers ending with .0 (like gp.Extent does)
    strExtent = (strExtent + " ").replace(".0 "," ").strip()

    return strExtent


def SetSnapRaster(SnapRasterDS,XReg=None,YReg=None,CellSize=None):
    """Sets the geoprocessor SnapRaster environment

    arguments

    SnapRasterDS - snap raster path (or raster layer)

    If the raster does not exist, a new one is created
    using the following parameters.

    XReg, YReg, CellSize - XY registration point and cell size

    These are ignored if SnapRasterDS exists.

    notes

    The current effective (snapped) extent is returned as a string.
    If the current extent isn't explicitly set (for example: "MAXOF",None)
    that's what is returned.

    This tool checks out a spatial analyst license
    if one is not already available.

    """
    global gp
    if not gp: gp = getGP()

    try:
        gp.CheckOutExtension("spatial")  # need a spatial analyst license
    except:
        raise Exception("Spatial Analyst license is not available")
    try:
        d = gp.Describe(SnapRasterDS)
        if d.DatasetType.find("Raster") == -1:
            raise Exception("%s is not a raster data set" % RasterDataset)
        SnapRasterDS = d.CatalogPath
        if XReg != None:
            GPMsg("w","Using existing Snap Raster %s" % SnapRasterDS)
    except Exception as xmsg:
        # Snap path does not exist, create one
        if not ( XReg and YReg and CellSize):
            raise Exception(
                  "XReg, YReg, and CellSize must be specified to create snap raster")
        # save environment
        tmpExtent = gp.Extent
        tmpCell = gp.CellSize
        gp.Extent = None
        gp.CellSize = CellSize
        gp.Extent = "%s %s %s %s" % \
          (XReg,YReg,XReg + CellSize * 2.1, YReg + CellSize * 2.1)
        gp.SnapRaster = None
        # Single Output Map Algebra will honor this extent always!
        gp.SingleOutputMapAlgebra("1.0",SnapRasterDS)
        gp.Extent = tmpExtent
        gp.CellSize = tmpCell

    # Set SnapRaster environment
    gp.SnapRaster = SnapRasterDS

    # Snap current extent to new Snap Raster
    try:
        gp.Extent = gp.Extent + " " + SnapRasterDS
    except:
        pass

    return gp.Extent


def DeleteList(FileList,Verbose=False):
    """Delete a ";"-delimited or python list of items

    (modified from ESRI HelperFunctions.py)

    example

      DeleteList("E:/work/temp1;e:/work/temp2")
      DeleteList(["file1","file2","file3"])

    """
    global gp
    if not gp: gp = getGP()

    if type(FileList) == type(""):
        FileList = FileList.split(";")

    for ff in FileList:
        if ff: # skip None or "" that may be in the list
            try:
                gp.Delete(f)
                if Verbose: GPMsg("Deleted \"%s\"" % ff)
            except:
                if Verbose: GPMsg("w","Could not delete \"%s\"" % ff)
                continue


def CalcField(in_table, field, expression, expression_type="PYTHON_9.3",
              code_block=None, ftype=None, flen=None):
    """CalculateField_management, improved

    improvements

    - Syntax defaults to PYTHON_9.3
    - If the field does not exist, it is added to the table first
    - If a field is added, Scale/Precision/Nullable/Required defaults are used

    arguments

    in_table (Table View / Raster Layer / Raster Catalog Layer / Mosaic Layer):

        The table containing the field that will be updated with the new calculation.

    field (Field):

        The field that will be populated with the calculation results.
        If the field does not exist, it will be added to the table first.

    expression (Python expression):

        Python expression to populate the selected rows.

    expression_type {String}:
        Specify the type of expression that will be used.
        Default is PYTHON_9.3 (best for 10.1 and later)

    code_block {String}:

        Allows for a block code to be entered for complex expressions.

    ftype {String}:

        Field type: ("LONG","DOUBLE",..) (See the tool help for Add Field)
        (Ignored if the field already exists.)

    flen {Long}:
        The length of the field being added.
        This option is only applicable on fields of type TEXT or BLOB.
        (Ignored if the field already exists.)

    """
    global gp
    if not gp: gp = getGP()

    if not gp.ListFields(in_table, field):
        gp.AddField_management(in_table, field,
                                  ftype, field_length=flen)
        #GPMsg("Added field {0} to {1}".format(field, in_table))

    gp.CalculateField_management(in_table, field, expression,
                                    expression_type, code_block)


def GetCount(input):
    """GetCount in one step (like 9.2 GetCount tool - except returns int)"""
    global gp
    if not gp: gp = getGP()

    Result = gp.GetCount(input)
    return int(Result.GetOutput(0))


def ListUnique(tbl, field, list_type="LIST"):
    """Create a list of unique values from a table/tableview.

    arguments

    tbl     Table or table view
    field   Field name
    list_type
      "LIST" - Python list of values
      "TEXT" - comma-separated list (i.e. for use in where expressions)
    """
    global gp
    if not gp: gp = getGP()
    try:
        try:
            # this will only work for 10.1
            field_values = [r[0] for r in
                            arcpy.da.SearchCursor(tbl, [field])]
        except:
            # ArcGIS 10.0
            Rows = gp.SearchCursor(tbl, "", "", field, field)
            Row = Rows.next()
            field_values = []
            while Row:
                field_values.append(Row.getValue(field))
                Row = Rows.next()
            del Row, Rows

        # unique-ize and sort the list
        field_values = sorted(set(field_values))

        # convert to text (comma-delimiated list) if requested
        if list_type.upper() == "TEXT":
            field_values = ",".join([repr(f) for f in field_values])
            # remove unicode prefix
            if field_values[0] == "u":
                field_values = field_values.replace("u'", "'", 1)
                field_values = field_values.replace(",u'", ",'")
        return field_values
    except:
        raise

def AdjustExtent(extent, amount=5, units="percent"):
    """Adjust an xy extent (default: +5 percent)

     extent - an extent object or string ("0 0 100 100")
     amount - amount to grow (+) or shrink (-) the extent
     units - "percent" or "mapunits"

    (modified from ESRI script HelperFunctions.py)

    NOTE: This has been superseded by the egis.GetExtent function.

    """
    # (modified from ESRI HelperFunctions.py, GetLargerExtent)
    # Adjust string Extent

    # AdjustExtent("1000 1000 2000 2000",5) == '950.0 950.0 2050.0 2050.0'
    # AdjustExtent("0 0 100 100",3,"map") ==  '-3.0 -3.0 103.0 103.0'
    # AdjustExtent("0 0 500 500",-5) == '25.0 25.0 475.0 475.0

    try:
        # convert extent to list of floats
        try:
            minX, minY, maxX, maxY = [ float(i) for i in strExtent.split() ]
        except:
            try:
                minX, minY, maxX, maxY = ext.Xmin, ext.YMin, ext.XMax, ext.Max
            except:
                raise Exception("Invalid extent: {}".format(extent))

        # if percent, (default) compute in map units
        if units == "percent":
            percent = float(amount) / 100.0
            amount = max(maxX - minX,maxY - minY) * percent

        # expand extent
        minX = minX - amount
        minY = minY - amount
        maxX = maxX + amount
        maxY = maxY + amount

        NewExtent = [minX, minY, maxX, maxY]
        NewExtent = [ str(i) for i in NewExtent ]
        return " ".join(NewExtent)
    except Exception as xmsg:
        GPMsg("AdjustExtent: Extent not adjusted\n"+\
              str(xmsg))


class MsgError(Exception):
    """Base class for errors in the egis module

    Example

      raise MsgError, "Houston, we have a problem...."

    """


def StringIsTrue(strText=""):
    """Converts value to boolean

      True only if the first character of string is:
        "t" ("true","totally","TRUTHINESS")
        or
        "y" ("yes", "yup", "Yowsa!")
        or argument is True or 1.
    """
    if str(strText).strip().lower()[0] in "1ty": return True
    else: return False


def TestValue(val, trueKey="true", falseKey="false", default=None):
    """Test a value and return a boolean result

    This function is useful to interpret booleans
    passed using keywords or other strings.
    Other data types are checked using the bool() function.

    arguments

    val       value to test
    trueKey   string keyword for True  (case-insensitive)
    falseKey  string keyword for False (case-insensitive)
    default   boolean to return if val is "#"

    example usage:

    overlap = sys.argv[2]
    overlap = TestValue(overlap, "OVERLAP", "NO_OVERLAP", True)
    """
    if isinstance(val, basestring):
        # convert string keywords to a boolean value
        val = val.strip().lower()
        if val in ["true", trueKey.lower(), "1", "-1"]:
            val = True
        elif val in ["false", falseKey.lower(), "", "0"]:
            val = False
        elif val == "#":
            if default != None:
                val = default
            else:
                raise Exception("No default specified for value '#'")
        else:
            raise Exception("Invalid keyword: '{0}'".format(val))
    # return boolean value
    return bool(val)


class GPModeThing:
    """A little Python class to keep track of print mode for GPMsg()

    Arguments

      mode -  "gp","print","both", or "off"
      logfile -  path to log file to echo messages

    See the help for the GPMsg() function for details.
      """
    def __init__(self):
        self.mode = "gp"
        self.logfile = None # default: no logfile

    def __call__(self, mode=None, logfile=None):
        #  "" strMode returns current value"""
        if mode:
            # check argument to make sure it is valid
            mode = mode.lower()
            if mode not in ["gp", "print", "both", "off"]:
                GPMsg("w", 'Valid values are: "gp","print","both", or "off"')
            else:
                self.mode = mode
        if logfile:
            if logfile.lower() == "close":
                self.logfile = None
            self.logfile = logfile
        return self.mode

# initialize it
GPMode = GPModeThing()


def GPMsg(sev=None, msg=None, msgType=None):
    """Send a message to the geoprocessor,"python-print", or both

    Geoprocessing messages displayed with methods like "gp.AddMessage"
    may be visible in ArcGIS or in tty python output, but when
    running the script within IDE environments like IDLE or Wing,
    invisible, or display garbled. This method allows you
    to specify to send a message to python-print, just the geoprocessor, or
    to both places. The syntax also allows for a little less typing
    than the gp messaging methods it calls.

    dependencies

    GetGP, GPModeThing

    arguments

    sev - severity code / message option

      "Message","Warning","Error","Return","Time","Gpmessages"
      (only the first letter is required)

     msg - text message to display

      A special syntax for msg is used to support gp.AddIDMessage().
      (Spaces are required between each argument!)

        ID <MessageID> {AddArgument1} {AddArgument2}

      For example, to do this:
        gp.AddIDMessage("Error", 12, outFeatureClass)
      You can use this syntax with GPMsg:
        GPMsg("Error","ID %s %s" % (12,outFeatureClass))
      If a message argument contains a space, you can use | to separate
        the second argument so it will parse correctly:
         GPMsg("Error","ID %s %s|%s" % (328,"Input dataset",outFeatureClass))
      (Please only use error numbers documented in the ArcGIS help!)

     msgType - Where to send the message. If this argument is given
        the destination will stay the same until GPMode is used or
        the argument is given again.

        "gp"      Geoprocessor (default) (arcpy or 9.x gp object)
        "print"   Python print
        "both"    Both places
        "off"     Nothing prints anywhere (use with care)
        None      Use current value of GPMode()

    examples

       GPMode("print")  # change default output to python-print
       GPMode(logfile="c:\\users\\jwpwell\\logs.txt") # echo messages to logfile
       GPMsg("This is a message") # default output, print ONLY
       x = egis.timeDiff() # reset timer
       GPMsg("t","The time is now","gp") # output to ArcGIS ONLY
       GPMsg() # print gp.AddMessages(0) GP messages to ARCGIS only
       GPMsg("w","ID 345","both") # use gp.AddIDMessage, output to both
       GPMode("off") # no messages printed
       # This code is used to pass along messages from a tool (skipping start and end times):
       arcpy.AddField_management("tbl", "FIELD1")
       for k in range(2, arcpy.GetMessageCount() - 1):  GPMsg("r", k)

       Output:

       This is a message
       10:40:05 2.34 The time is now
       Executing: CopyFeatures poly_inout E:\work\poly_inout_CopyFeatures.shp # 0 0 0
       Start Time: Wed Apr 07 11:11:58 2010
       Executed (CopyFeatures) successfully.
       End Time: Wed Apr 07 11:11:58 2010 (Elapsed Time: 0.00 seconds)
       WARNING 000345: Input must be a workspace
       Adding FIELD1 to tbl...

    """
    global gp
    if not gp: gp = getGP()

    # support shorthand usage: GPMsg("message") and GPMsg()
    if sev != None and msg == None:
        # GPMsg("message") ->  GPMsg("","message")
        sev,msg = None,sev
    elif sev == None and msg == None:
        # GPMsg() -> GPMsg("Message",gp.GetMessages(0))
        sev,msg  = "g",None

    if not msgType:
        msgType = GPMode()  # use current value of GPMode
    else:
        msgType = GPMode(msgType) # set GPMode (and remember for next GPMsg)

    if msgType == "off":
        # Do not print anything! (like AML "&messages &off &all")
        return

    # decode severity to a code 0 thru 5:
    # sev  isev description
    #        0  message
    # "w"    1  warning
    # "e"    2  error
    # "r"    3  return message  (next argument is message #)
    # "t"    4  message with clock time
    # "g"    5  return gp.GetMessages(0)
    dictSev = {"w":1, "e":2, "r":3, "t":4,"g":5 }
    try:
        sev = str(sev).lower()[:1] # "Warning" -> "w"
        isev = dictSev[sev]
    except:
        isev = 0

    # support gp.AddIDMessage
    IDMessage = False # assume this isn't an id message
    try:
        # Usage: "ID <msgID> {|} {arg1} {|} {arg2}"
        lstMsg = msg.split()
        if lstMsg[0].lower() != "id": raise
        try:
            MessageID = int(lstMsg[1])
            if MessageID <= 0 or MessageID > 99999: raise
        except:
            GPMsg("w","GPMsg: Invalid message ID: %s" % MessageID)
            raise
        xmsg = " ".join(lstMsg[2:])
        if xmsg.find("|") > -1:
            lstMsg = xmsg.split("|")
        else:
            lstMsg = xmsg.split()

        IDMessage = True
        # capture AddArgument1, AddArgument2
        try:
            IDArg1, IDArg2 = "", ""
            IDArg1 = lstMsg[0]
            IDArg2 = lstMsg[1]
        except:
            pass
    except:
        pass

    # if time asked for, calc time tag
    if isev == 4:
        curTime = time.strftime("%H:%M:%S", time.localtime())
        elpTime = timeDiff()
        msg_time = "%s %05.2f" % (curTime, elpTime)
    # if asking for geoprocessing message only, capture it
    if isev == 5:
        msg = gp.GetMessages(0)

    # send our message

    if msgType.lower() in ["gp","both"]:
        # send message to geoprocessor
        if isev == 0:
            gp.AddMessage(msg)
        elif isev == 1:
            if not IDMessage:
                gp.AddWarning(msg)
            else:
                gp.AddIDMessage("Warning",MessageID,IDArg1,IDArg2)
        elif isev == 2:
            if not IDMessage:
                gp.AddError(msg)
            else:
                gp.AddIDMessage("Error",MessageID,IDArg1,IDArg2)
        elif isev == 3:
            gp.AddReturnMessage(int(msg))
        elif isev == 4:
            gp.AddMessage("%s %s" % (msg_time, msg))
        elif isev == 5:
            gp.AddMessage(msg)
    if msgType.lower() in ["print","both"]:
        # python-print messages
        SevLabel = ["","WARNING","ERROR"]
        if isev == 0:
            print(msg)
        elif isev in [1,2]:
            if not IDMessage:
                print("%s: %s" % (SevLabel[isev],msg))
            else:
                print("%s %06d: %s" % (SevLabel[isev],MessageID,lstMsg))
        elif isev == 3:
            print(gp.GetMessage(int(msg)))
        elif isev == 4:
            print("%s %s" % (msg_time, msg))
        elif isev == 5:
            if msg:
                print(msg)

    # write message to log if set
    if GPMode.logfile:
        logf = open(GPMode.logfile, "a")
        SevLabel = ["","WARNING","ERROR"]
        if isev == 0:
            logf.write("%s" % msg)
        elif isev in [1,2]:
            if not IDMessage:
                logf.write("%s: %s" % (SevLabel[isev], msg))
            else:
                logf.write("%s %06d: %s" % (SevLabel[isev],MessageID,lstMsg))
        elif isev == 3:
            logf.write(gp.GetMessage(int(msg)))
        elif isev == 4:
            logf.write("%s %s" % (msg_time, msg))
        elif isev == 5:
            if msg:
                logf.write(msg)
        logf.write("\n")
        logf.close()


def TraceInfo():
    """Returns traceback information for easy reporting.
       Modified from Dale Honeycutt's (ESRI) blog post:
       "Error handling in Python script tools":
       http://blogs.esri.com/esri/arcgis/2008/12/01/tips-and-tricks-error-handling-in-python-script-tools/
       Curtis Price - cprice@usgs.gov - 2009-05-06

      Here's an example of how to use it:

         except:
           line, file, err = egis.TraceInfo()
           gp.AddError("Python error on %s of %s\n" % (line,file,err)
           # or
           arcpy.AddError("Python error on %s of %s\n" % (line,file,err)
    """
    import sys, os, traceback
    try:
        # get traceback info
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]  # script name + line number
        # Get line number
        ErrLine = tbinfo.split(", ")[1]
        # Get error message
        ErrMsg = traceback.format_exc().splitlines()[-1]
    except:
        # just in case *this*  thing fails!
        Here = os.path.realpath(__file__)
        Msg = traceback.format_exc().splitlines()[-1]
        Msg = "TraceInfo failed\n%s\n(%s)" % (Here, Msg)
        raise Exception(Msg)
    return ErrLine, sys.argv[0], ErrMsg


def CheckOutExtension(strExtName):
    """Checks out an ArcGIS extension.

    Like gp/arcpy.CheckOutExtension(), except returns the string
    "NotInstalled" if the requested extension is not installed.

    arguments

      strExtName - Case-insensitive extension name. Valid names include:
        "3d","datainteroperability","geostats","mpsatlas","network",
        "schematics","spatial","streetmap","survey","tracking"

    examples

      strStat = egis.CheckOutExtension("spatial")
    """
    global gp
    if not gp: gp = getGP()

    import os
    try:
        # Set up dictionary of "tag paths" to identify whether
        # each extension by whether a certain file exists

        TagPath = {\
            "spatial"             : "Bin/GridCore.dll",\
            "3d"                  : "Bin/3DAnalystUI.dll",\
            "streetmap"           : "Locators/US Single House.avs",\
            "mpsatlas"            : "Solutions/MPSAtlas",\
            "network"             : "NetworkAnalyst",\
            "schematics"          : "Schematics",\
            "survey"              : "SurveyAnalyst",\
            "datainteroperability": "Data Interoperability Extension",\
            "geostats"            : "Bin/GATools.dll",\
            "tracking"            : "Bin/TrackingCore.dll" \
        }
        lstTag = TagPath.keys()
        lstTag.sort()
        strTags = ','.join(["\"%s\"" % k for k in lstTag])

        # is extension installed?
        strExt = strExtName.lower()
        if strExt in TagPath:
            ARCGISHOME = gp.GetInstallInfo("desktop")["InstallDir"]
            strFlag = os.path.join(ARCGISHOME,TagPath[strExt])
        else:
            strMsg = "ArcGIS extension \"%s\" not supported.\n" % strExt
            strMsg += "Valid keys: " + strTags
            raise Exception(strMsg)
        if not os.path.exists(strFlag): return "NotInstalled"

        # OK, it's installed, go ahead and try to check it out
        return gp.CheckOutExtension(strExt)

    except Exception as xmsg:
        raise Exception(str(xmsg))


def ListPaths(inDir=".",strWild=".*",FoldersOnly=False,FullRegExp=False):
    """Return a list of file or folder paths using Python os.walk function

    This function is used by the Lister script tool in the EGIS toolbox.
    """
    import os, re
    lstOut = []
    if not FullRegExp:
        # convert simple regexp to full regexp
        strWild0 = strWild
        strWild = strWild.replace(r".",r"\.")
        strWild = strWild.replace(r"*",r".*")
        strWild += r"$"
    GPMsg("Matching regular expression: \"%s\" ..." % strWild)
    regMatch = re.compile(strWild,re.IGNORECASE)
    try:
        for root, dirs, files in os.walk(inDir):
            if regMatch.match(root) != None:
                lstOut.append(os.path.realpath(root))
            if not FoldersOnly:
                for strFileName in files:
                    # find a match?
                    if regMatch.match(strFileName) != None:
                        lstOut.append(\
                            os.path.realpath(os.path.join(root,strFileName)))
        return lstOut
    except Exception as ErrorDesc:
        raise Exception("Search failed for \"%s\"\n" % strWild + str(ErrorDesc))


def SysCommands(commandString=None,RunFolder=".",WarnMsg=False):
    """Runs a command shell

    This tool sends commands to a command shell and prints dialog

    arguments

      commandString - text string of tty to feed to command process

      RunFolder - Folder in which to start up command process

      WarnMsg - boolean - display stdout messages as warnings

    examples

     import egis
     Sev = egis.SysCommands("dir","e:/work")

    (Note: stderr is captured as GP error messages [cool, huh?])
    """
    global gp
    if not gp: gp = getGP()

    try:
        # Command string
        if not commandString:
            raise MsgError("No command provided")

        # does the workspace folder exist?
        if RunFolder == "#" or RunFolder == "":
            RunFolder = os.environ["TEMP"]
        if not os.path.exists(RunFolder):
            GPMsg("w","\"%s\" not found" % RunFolder)
            GPMsg("w","ID 535")
            RunFolder = os.environ["TEMP"]
            GPMsg("w","Using TEMP folder: %s" % RunFolder)

        # Submit the command to the shell
        lstCmd = commandString.split("\n")

        if len(lstCmd) > 1:
            # prepend multi line commands with a cmd /c
            lstCmd = ['cmd /c'] + lstCmd
        # else:
            # Single commands are run with no dialog
            # GPMsg(commandString)

        # Go to folder, start shell
        strHere = os.path.realpath(os.curdir) # save it
        os.chdir(RunFolder)
        fi,fo,fe = os.popen3(lstCmd[0],'t')

        # for multi-line command session, send each command one at a time
        if len(lstCmd) > 1:
            for Cmd in lstCmd[1:]:
                fi.write(Cmd + "\n")
            # add command to close the shell
            fi.write("exit\n")

        # Read the session dialog and pass it along to the GP.
        # The data are not blocked so they will be reported
        # as they come back to the GP message stream.
        strMsg = fo.readline()
        Sev = 0 # no errors so far!

        if WarnMsg: MsgType = "w"
        else: MsgType = ""
        # read std output
        while strMsg != "":
            GPMsg(MsgType,strMsg.strip())
            strMsg = fo.readline() # Get next line of output

        # check for std error messages
        strErr = fe.readline()
        if strErr != "": Sev = 2 # we have an error of some kind
        while strErr != "":
            GPMsg("E",strErr.strip())
            strErr = fe.readline() # Get next line of error

        # clean up
        fi.close()
        fo.close()
        fe.close()
        del fi,fo,fe
        os.chdir(strHere)
        return Sev

    except Exception as ErrorDesc:
        # Return error messages
        GPMsg("e",str(ErrorDesc))

        try:
            os.chdir(strHere)  # go back to this folder once shell started
        except:
            pass
        # return error severity, in case there is any doubt
        return 2


def ArcCommands(commandString, ArcWorkspace=None, strEcho="off"):
    """Run an ArcInfo Workstation session

    arguments

      commandString - ArcInfo Workstation commands

      ArcWorkspace - Folder to start session

      strEcho -"&echo" setting: ("off" [default], "on", "brief")

    examples

      Sev = ArcCommands("clean mycov # # # 1.0;labelerrors mycov",r"e:\work")

    notes

    Take care that the dialog specified in commandString avoids stopping
    the session inside a dialog - if so, the only way to kill the workstation
    session is by killing the arc.exe process.

    The current geoprocessing environment coverage settings are used.

    """
    global gp
    if not gp: gp = getGP()

    try:
        # Is Workstation installed?
        try:
            ARCHOME = os.environ["ARCHOME"]
        except:
            raise MsgError("ArcInfo Workstation is not installed.")

        # Command string
        if commandString == "#" or commandString.strip() == "":
            raise MsgError("No command provided")

        # check command lines for long length (ArcInfo limit)
        lstCommands = commandString.split("\n")
        k = 0
        for strCmd in lstCommands:
            if len(strCmd) > 320:
                strMsg = "Input line %s too long (%s/320 max)\n\n%s\n\n" \
                       % (str(k),str(len(strCmd)),strCmd)
                raise MsgError(strMsg)
            k += 1

        # does the workspace folder exist?
        if not ArcWorkspace: ArcWorkspace = gp.Workspace
        if not ArcWorkspace: ArcWorkspace = os.environ["TEMP"]

        if not os.path.exists(ArcWorkspace):
            GPMsg("w","Workspace \"%s\" not found" % ArcWorkspace)
            ArcWorkspace = os.environ["TEMP"]
            GPMsg("w","Using %s" % ArcWorkspace)
        elif ArcWorkspace.find(".mdb") > 0 or ArcWorkspace.find(".gdb") > 0:
            GPMsg("w","Workspace must be a filesystem folder (not a geodatabase)")
            ArcWorkspace = os.path.dirname(ArcWorkspace)
            GPMsg("w","Using %s" % ArcWorkspace)

        ArcWorkspace = os.path.dirname(GetShortName(ArcWorkspace + "/x"))

        # ArcInfo Workstation &echo
        strEcho = strEcho.upper()
        if strEcho.find("OFF")>-1: strEcho = "off"
        elif strEcho.find("BRI")>-1: strEcho = "brief"
        elif strEcho.find("ON")>-1: strEcho = "on"
        else: strEcho = "off"

        # GP Coverage setting defaults
        if not gp.NewPrecision: gp.NewPrecision = "DOUBLE"
        if not gp.DerivedPrecision: gp.DerivedPrecision = "HIGHEST"
        if not gp.ProjectCompare: gp.ProjectCompare = "NONE"
        strGPEnv = "precision %s %s;projectcompare %s" % \
                 (gp.NewPrecision.lower(), gp.DerivedPrecision.lower(),\
                  gp.ProjectCompare.lower())

        # Go to folder
        strHere = os.path.realpath(os.curdir)

        os.chdir(ArcWorkspace) # run arc shell from this workspace

        # command to start ArcInfo Workstation shell
        ArcCmd = os.path.realpath(ARCHOME + "/bin/arc.exe")  + "\n"

        # Start an interactive ArcInfo session and feed it the command line
        # we need the 'b' argument to properly get the newlines
        fi,fo,fe = os.popen3(ArcCmd)

        # Send interactive commands to the Arc prompt

        # set up session
        # Geoprocessing Coverage settings
        strCommand = strGPEnv + "\n"
        # ignore errors (we will trap them below)
        strCommand += "&severity &warning &ignore;&severity &error &ignore"
        # set &echo environment if the user asked for it
        if strEcho != "off":
            strCommand += ";&echo &" + strEcho.lower()
        fi.write("%s\n" % strCommand)

        # Send the command the user asked for
        fi.write("%s\n" % commandString)

        # report severity in a message so we can capture the error status
        # and pass it along to the geoprocessor later
        fi.write("&type ***AML Severity: %aml$sev%, " +\
                 "%aml$message%\n")

        # Send the "quit" command to end the session, but with
        # some insurance... stack in extra newlines and quit commands
        # just in case so we do our best to end the shell
        fi.write("quit\nquit\n" + "\n" * 10 + "quit\nquit\n")


        # We've sent all we can send to the arc process...

        # Now, read the session dialog and pass it along to the GP.
        # The data are not blocked so they will be reported
        # as it happens to the GP message stream.

        strMsg = str(fo.readline())

        # print version startup lines and brief license notice
        for line in range(3):
            GPMsg(strMsg.strip())
            strMsg = fo.readline()
        # skip startup verbiage until arc prompt
        while not strMsg.startswith("Arc:"):
            strMsg = fo.readline()

        # initialize maxSeverity found
        maxSeverity = 0

        # parse output

        while strMsg != "":
            # strip newlines
            strMsg = strMsg.strip()

            # Check output, and echo results
            Severity = ArcMsgSeverity(strMsg)

            # replace double \\ with \n:
            strMsg = strMsg.replace(r"\\","\n: ")

            # modify severity report string
            if strMsg.find("***AML") > -1:
                strMsg = strMsg.replace("***AML Severity: 1, ","")
                strMsg = strMsg.replace("***AML Severity: 2, ","")
            if Severity == 0:
                pass # in the skip list
            elif Severity == 1:
                GPMsg("w",strMsg)
            elif Severity == 2:
                GPMsg("e",strMsg)
            else: # ok to print
                GPMsg("",strMsg)

            # record max severity code encountered
            maxSeverity = max(maxSeverity,Severity)

            strMsg = fo.readline() # Get next ArcInfo session message


        # check for std error messages
        strErr = fe.readline()
        if strErr != "": maxSeverity = 2 # we have a system error
        while strErr != "":
            GPMsg("E",strErr[:-1])
            strErr = fe.readline() # Get next line of error

            # clean up
        fi.close()
        fo.close()
        del fi, fo, fe
        os.chdir(strHere)  # go back to this folder once shell started

        return maxSeverity
    except MsgError as xmsg:
        # this is an error to pass on to the user
        raise MsgError(str(xmsg))
    except arcgisscripting.ExecuteAbort(xmsg):
        # user hit cancel button
        GPMsg("e",str(xmsg) + \
              "\nYou may have to manually kill the process \"arc.exe\"" + \
              "from Windows Task Manager.")
    except Exception:
        # Python errors
        line,file,err = TraceInfo()
        raise MsgError("Python error on %s of %s:\n%s" % (line,file,err))
    finally:
        try:
            os.chdir(strHere)  # go back to this folder once shell started
        except:
            pass


def ArcMsgSeverity(strMessage):
    """Scans a string for ArcInfo messages, returning a severity code
    """
    # this routine checks input/output messages from the arcinfo session
    # dialog for strings that tell us error messages. This function
    # returns a severity code:
    #  0 = no warning or error
    #  1 = this is a warning
    #  2 = this is an error
    # -1 = no matches found

    # input/output strings to skip make the output cleaner
    SkipList = ["Arc: &severity &warning &ignore;&severity &error &ignore",
                "&type ***AML Severity:",
                "***AML Severity: 0",]
    # Warning and Error strings
    # "FATAL" and "Submitting command to operating system" message are flagged
    # to generate warning and errors (in AML they don't)
    WarnList = [
        "Submitting command to Operating System",
        "***AML Severity: 1,",
        "AML WARNING"
    ]
    ErrList = [
        "FATAL",
        "***AML Severity: 2,",
        "AML ERROR"
    ]

    Sev = 0
    for List in (SkipList, WarnList, ErrList):
        for strMsg in List:
            if strMessage.find(strMsg) >= 0:
                return Sev
        Sev = Sev + 1
    return -1  #  no matches found, return -1


def GetShortName(strFullPath):
    """Convert the directory portion of a path string to short filename style.

    Note that it is still possible for the name to be too long or to
    have some character in it that is not supported by Workstation.

    arguments

      strFullPath - A full pathname

    example

    >>> egis.getShortName("c:/Documents and settings/cprice/Desktop")
    'c:\\DOCUME~1\\cprice\\Desktop'

    See the ArcGIS help for a discussion of paths in ArcGIS.

    This function was modified from the Arc 10.0 script tool
    Desktop10.0\ArcToolbox\Scripts\importe00.py
    """
    import subprocess

    # TODO: If Unix, just return the input, with a warning
    # if the path has spaces in it.

    # if this is an existing folder path with no trailing os.sep, add it
    if os.path.isdir(strFullPath):
        if strFullPath[-1] not in ["\\","/"]:
            strFullPath += os.sep

    (long_dir, long_name) = os.path.split(strFullPath)

    for_cmd = 'for %I in ("' + long_dir + '") do echo %~sI'
    p = subprocess.Popen(for_cmd, shell=True, stdout=subprocess.PIPE).stdout
    short_dir = p.readlines()[-1] # last line from for command
    if p.close():
        #This means there was an error calling shell command "for"
        #
        #If the incoming full path is short enough and has no spaces,
        # and the incoming name is short enough, then we can just use
        # it, otherwise, we need to throw an error. (max length = 122)
        if ( (len(strFullPath)<123)
           and (strFullPath.find(" ") == -1)
           and (len(strFullPath)<=13) ):
            return (strFullPath)
        else:
            # 001103 : Error converting directory to short name format.
            arcpy.AddIDMessage("Error", 1103)
            raise Exception

    #Strip whitespace off the end
    short_dir = short_dir.rstrip()

    #Add the unshortened file portion back onto the now-shortened path
    short_path = os.path.join(short_dir, long_name)

    return (short_path)

def pprint_fields(table, wild=None, detail=False):
    """ pretty print table's fields and their properties """
    def _print(l):
        GPMsg("".join(["{:>12}".format(i) for i in l]))

    atts =  ['name', 'aliasName', 'type', 'baseName']
    if detail:
        atts += ['domain', 'editable', 'isNullable', 'length', 'precision',
            'required', 'scale']
    d = arcpy.Describe(table)
    GPMsg("{} {}".format(d.dataType, d.catalogPath))
    _print(atts)
    for f in arcpy.ListFields(table, wild):
        _print(["{:>13}".format(getattr(f, i)) for i in atts])

TEMPLATE = """
# Template script for a an ArcPy script tool using egis.py
# Shows how to use the egis module to trap errors,
# manage messages, and handle scratch files.
#
# Use this script as a starting point for your own scripting.
#
# More geoprocessing script template ideas:
# http://blogs.esri.com/esri/arcgis/2011/08/04/pythontemplate/
# http://blogs.esri.com/esri/arcgis/2008/12/01/tips-and-tricks-error-handling-in-python-script-tools/

# History:
# Clarence King, cking@usgs.gov, 1879-03-03, original coding
# Curtis Price, cprice@usgs.gov, 2011-10-12 arcpy compatibility; fixes
# Curtis Price, 2011-11-15 updated exception handling
# Curtis Price, 2014-05-01 updated exception handling to work with python 3x
# Curtis Price, 2014-11-14 tool encapsulated in function

# Python 3.x compatibility
from __future__ import print_function, unicode_literals, absolute_import

import sys
import os
import traceback
import arcpy
## arcpy.CheckOutExtension("spatial") # spatial analyst license
## from arcpy.sa import * # for ArcPy map algebra
import egis
# from/import commonly used egis methods
from egis import GPMsg, MsgError

def test_tool(arg=None):
    try:
        # egis.ScratchName method gets valid scratch names
        tmpGrid = egis.ScratchName("","","grid")

        # Three different kinds of errors
        # (uncomment one at a time to see how they get handled)

        # 1. Script-raised error - send a message to user
        #raise MsgError, "Your data have issues"
        # 2. runtime Python error, divide by zero
        #print 3 / 0
        # 3. ESRI Geoprocessing error (grid doesn't exist)
        #Result = arcpy.GetCount_management(tmpGrid)

        GPMsg("If this you see this, there were no errors!")
        GPMsg("w", "This is warning, just to keep you on your toes.")
        # ArcGIS message w/ args
        GPMsg("w", "ID 732 %s|%s" % ("tmpGrid",tmpGrid))
        # script timer
        GPMsg("t", "Now is the time for all good persons...")

    except MsgError as xmsg:
        GPMsg("e", str(xmsg))
    except arcpy.ExecuteError:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        GPMsg("e", tbinfo.strip())
        for i in range(arcpy.GetMessageCount()):
            GPMsg("r", i) # use AddReturnMessage()
    except Exception as xmsg:
        # generic python errors
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        GPMsg("e", tbinfo.strip())
        GPMsg("e", str(xmsg))
    finally:
        # Clean up here (delete cursors, temp files)
        pass

if __name__ == "__main__":
    # Script tool interface:
    # Get zero or more parameters as text, call tool
    argv = tuple(arcpy.GetParameterAsText(i)
        for i in range(arcpy.GetArgumentCount()))
    test_tool(*argv)
"""

TEMPLATE9 = """
# Template script for a GP script tool using egis.py
# Shows how to use the egis module to trap errors,
# manage messages, and handle scratch files.
#
#------------------------------------------------------------------
# ArcGIS 9.x version (arcgisscripting)
#------------------------------------------------------------------
#
# Use this script as a starting point for your own scripting.
#
# More geoprocessing script template ideas:
# http://blogs.esri.com/esri/arcgis/2011/08/04/pythontemplate/
# http://blogs.esri.com/esri/arcgis/2008/12/01/tips-and-tricks-error-handling-in-python-script-tools/
#
# History:
# Clarence King, crking@usgs.gov, 1879-03-03, original coding
#                http://en.wikipedia.org/wiki/Clarence_King

import sys
import os
import traceback

import egis
# from/import commonly used egis methods
from egis import GPMsg, MsgError

gp = egis.getGP(9.3,"spatial") # create gp; check out spatial analyst

try:

    # get arguments
    arg1 = gp.GetParameterAsText(0)

    # egis.ScratchName method gets valid scratch names
    tmpGrid = egis.ScratchName("","","grid")

    # Three different kinds of errors
    # (uncomment one at a time to see how they get handled)

    # 1. Script-raised error - send a message to user
    #raise MsgError, "Your data have issues"
    # 2. runtime Python error, divide by zero
    #print 3 / 0
    # 3. ESRI Geoprocessing error (grid doesn't exist)
    #Result = arcpy.GetCount_management(tmpGrid)

    GPMsg("If this you see this, there were no errors!")
    GPMsg("w", "This is warning, just to keep you on your toes.")
    GPMsg("w", "ID 732 %s|%s" % ("tmpGrid",tmpGrid)) # ArcGIS message w/ args
    GPMsg("t", "Now is the time for all good persons...") # script timer

# handle errors and report

except MsgError, xmsg:
    GPMsg("e", str(xmsg))
except arcgisscripting.ExecuteError:
    tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
    GPMsg("e", tbinfo.strip())
    for i in range(gp.MessageCount):
        GPMsg("r", i) # use gp.AddReturnMessage
except Exception,  msg:
    # generic python error
    tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
    GPMsg("e", tbinfo.strip())
    GPMsg("e", str(msg))
finally:
    # Clean up here (delete cursors, temp files)
    pass

"""
