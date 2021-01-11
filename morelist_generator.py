#!/usr/bin/env python
#####################
#
# list_generator
#
# Generates lists of production or
# processed data (full RAT format or
# ntuples) for download.
#
# Author: Matt Mottram
#         <m.mottram@sussex.ac.uk>
#
# Revision: Michal Rigan <m.rigan@sussex.ac.uk> 11/03/2019: Fix for `None` in the surl list for production_list generation with exactly file.
#
#####################

import os
import sys
import database
import grid
import utilities
try:
    import json
except:
    raise ImportError("json not available; python version needs to be 2.7+")
    

# Copy type will be set to SURL if lcg-utils or gfal are not available
copy_type = "GUID"
server_list = []
# The first allowed type will ALWAYS be the default
#allowed_file_types = ["ntuple", "ratds", "soc"]
allowed_file_types = ["ntuple", "ratds", "soc", "ratds_blind", "ratds_unblind"]

def get_module_list_version(data_type, rat_version=None):
    """Get a list of available modules.
    """
    rows = database.view("_design/proddata/_view/data_by_mod_pass_sr",
                         reduce=True, group=True, group_level=3,
                         startkey=[data_type], endkey=[data_type,{}])
    module_list = []
    for row in rows:
        if rat_version and row["key"][2] != unicode(rat_version):
            continue
        module_list.append(row["key"][1])
    return list(set(module_list))


def get_module_list_label(data_type, label=None):
    """Get a list of available modules.                                                                                                            
    """
    rows = database.view("_design/proddata/_view/data_by_label",
                         reduce=True, group=True, group_level=3,
                         startkey=[data_type,label], endkey=[data_type,label,{}])
    module_list = []
    for row in rows:
        if label and row["key"][1]!=unicode(label):
            continue
        module_list.append(row["key"][2])
    return list(set(module_list))

def rows_to_list(rows, run_list, all_passes, ratdb_tag):
    ''' 
    Takes a couchdb query result (rows) and turns it into names, sizes, guids, adlers.
    run_list specifies specific runs to consider (as a list, if not None)
    all_passes can either be boolean (True -> return all passes, False -> only newest)
      OR all_passes can be an integer specifying a specific pass to grab
    If run_list and all_passes are both None everything in rows will be considered.
    Finally, if ratdb_tag is not None, only files matching the specific database 
      tag will be returned (may be a string or list of strings), otherwise returns 
      all files satisfying previous conditions.
    '''
    if isinstance(ratdb_tag,str):
        ratdb_tag = [ratdb_tag]
    names = []
    sizes = []
    guids = []
    adlers = []
    pass_map = {}
    for row in rows:
        run = row["key"][-3]
        if run_list is not None and run not in run_list:
                continue
        
        # Get the details of this run
        name = row['value'][0]
        size = row['value'][1]
        if copy_type=="GUID":
            guid = row["value"][2]
            if not guid or guid=="":
                print "Cannot get file %s, no GUID" % (row["value"][0])
                continue
        else:
           # Fix for 'None' in surl list
            surl_list = row["value"][4]
            if (len(surl_list)>1):
                for lsurl in surl_list:
                        if lsurl is None:
                                surl_list.remove(None)
                                print "Removed `None` from the surl list"
            guid = grid.get_closest_copy(server_list, surl_list)
            if len(guid)==0 or guid=="":
                print "Cannot get file %s, no SURL" % row["value"][0]
                continue
        adler = row["value"][5]
        
        # Filter out old passes if desired
        if type(all_passes) is bool and not all_passes:
            passno = row['key'][-1]
            subrun = int(row['key'][-2])
            if (run,subrun) in pass_map:
                lastpassno,idx = pass_map.pop((run,subrun))
                if lastpassno < passno:
                    names[idx] = name
                    sizes[idx] = size
                    guids[idx] = guid
                    adlers[idx] = adler
                    pass_map[(run,subrun)] = (passno,idx)
                else:
                    pass_map[(run,subrun)] = (lastpassno,idx)
                continue #either replaced the older pass or this one was older
            else:
                pass_map[(run,subrun)] = (passno,len(names))
        elif type(all_passes) is int:
            passno = row['key'][-1]
            if passno != all_passes:
                continue #don't add a pass that is not the pass we want
        
        if type(ratdb_tag) is list:
            datadoc = database.table(row['id'])
            if 'ratdb_tag' not in datadoc:
                continue #has no database tag
            if passdoc['ratdb_tag'] not in ratdb_tag:
                continue #tag is not the tag we want
        
        names.append(name)
        sizes.append(size)
        guids.append(guid)
        adlers.append(adler)
    return names, sizes, guids, adlers

def generate_list_from_exactly(data_type, exactly):
    '''Gets a list of data from a list of module,run,pass (exactly).'''    
    names = []
    sizes = []
    guids = []
    adlers = []
    for module,run,passno in exactly:
        startkey = [data_type, module, run]
        endkey = [data_type, module, run, {}]
        rows = database.view("_design/proddata/_view/data_by_mod_run",
                       startkey=startkey, endkey=endkey, reduce=False)
        n, s, g, a = rows_to_list(rows,None,passno,None)
        names += n
        sizes += s
        guids += g
        adlers += a
    return names, sizes, guids, adlers

def get_list_by_module(data_type, module, run_range, run_list=None, all_passes=True, ratdb_tag=None):
    """Get a list of files for the settings given.
    """
    if run_range:        
        startkey = [data_type, module, run_range[0]]
        endkey = [data_type, module, run_range[1]] 
    else:
        startkey = [data_type, module]
        endkey = [data_type, module, {}]
    rows = database.view("_design/proddata/_view/data_by_mod_run",
                   startkey=startkey, endkey=endkey, reduce=False)

    return rows_to_list(rows,run_list,all_passes,ratdb_tag)


def get_list_by_module_version(data_type, module, rat_version, run_range, run_list=None, all_passes=True, ratdb_tag=None):
    """Get a list of files for the settings given.
    """
    if run_range:        
        startkey = [data_type, module, rat_version, run_range[0]]
        endkey = [data_type, module, rat_version, run_range[1]] 
    else:
        startkey = [data_type, module, rat_version]
        endkey = [data_type, module, rat_version, {}]
    rows = database.view("_design/proddata/_view/data_by_mod_ratv_run",
                   startkey=startkey, endkey=endkey, reduce=False)

    return rows_to_list(rows,run_list,all_passes,ratdb_tag)


def get_list_by_label_module(data_type, label, module, run_range, run_list=None, all_passes=True, ratdb_tag=None):
    """Get a list of files for the settings given

    If no module give, return all modules.
    """
    if run_range:        
        startkey = [data_type, label, module, run_range[0]]
        endkey = [data_type, label, module, run_range[1]] 
    else:
        startkey = [data_type, label, module]
        endkey = [data_type, label, module, {}]  
    rows = database.view("_design/proddata/_view/data_by_label_mod_run",
                   startkey=startkey, endkey=endkey, reduce=False)

    return rows_to_list(rows,run_list,all_passes,ratdb_tag)


def get_labels(data_type, module=None):
    """Get a list of all labels available

    Specify a module to see versions associated with that module.
    """
    labels = set()
    if module is None:
        # List any labels
        rows = database.view("_design/proddata/_view/data_by_label",
                             reduce=True, group_level=2,
                             startkey=[data_type], endkey=[data_type, {}])
        for row in rows:
            lb = row["key"][1]
            if type(lb)==list:
                lb = lb[0]
            labels.add(lb)
    elif type(module)==list and len(module)==1:
        # List labels only for a given module
        rows = database.view("_design/proddata/_view/data_by_label",
                             reduce=True, group_level=3,
                             startkey=[data_type], endkey=[data_type, {}])
        for row in rows:
            if str(row["key"][2])==module[0]:
                lb = row["key"][1]
                if type(lb)==list:
                    lb = lb[0]
                labels.add(lb)
    else:
        raise Exception("Must only provide one module name.")
    return labels


def get_versions(data_type, module=None):
    """Get a list of all versions available.

    Specify a module to see versions associated with that module.
    """
    versions = set()
    if module is None:
        # List all RAT versions
        rows = database.view("_design/proddata/_view/data_by_mod_pass_sr",
                             reduce=True, group_level=3,
                             startkey=[data_type], endkey=[data_type, {}])
        for row in rows:
            versions.add(row["key"][2])
    elif type(module)==list and len(module)==1:
        # List only versions for a given module
        rows = database.view("_design/proddata/_view/%s_by_mod_pass_sr",
                             reduce=True, group_level=3,
                             startkey=[data_type], endkey=[data_type, {}])
        for row in rows:
            if str(row["key"][1])==module[0]:
                versions.add(row["key"][2])
    else:
        raise Exception("Must only provide one module name.")
    return versions


def get_modules(data_type, label=None, version=None):
    """Get a list of all modules available

    Only one of label or version can be specified (or neither).
    """
    if label is not None and version is not None:
        raise Exception("May only specify label OR version")
    modules = set()
    if label is not None:
        # Only modules for a given label
        rows = database.view("_design/proddata/_view/data_by_label" % (data_type),
                             reduce=True, group_level=3, startkey=[data_type, label],
                             endkey=[data_type, label, {}])
        for row in rows:
            modules.add(row["key"][2])
    elif version is not None:
        # Only modules for a given RAT version
        rows = database.view("_design/proddata/_view/data_by_mod_pass_sr",
                             reduce=True, group_level=3, startkey=[data_type],
                             endkey=[data_type, {}])
        for row in rows:
            if str(row["key"][2])==version:
                modules.add(row["key"][1])
    else:
        # Any and all modules
        rows = database.view("_design/proddata/_view/data_by_mod_pass_sr",
                             reduce=True, group_level=2,
                             startkey=[data_type], endkey=[data_type, {}])
        for row in rows:
            modules.add(row["key"][1])
    return modules


def get_response(host, url, username=None, password=None):
    headers = {}
    if username is not None and password is not None:
        auth_string = base64.encodestring('%s:%s' % (username, password))[:-1]
        headers['Authorization'] = 'Basic %s' % auth_string
    connection = httplib.HTTPConnection(host, port=5984)
    try:
        connection.request('GET', url, headers=headers)
        response = connection.getresponse()
    except httplib.HTTPException as e:
        sys.stderr.write('Error accessing the requested db query: %s' % str(e))
        sys.exit(20)
    return response.read()


def generate_list_by_label(label, data_type, module, run_range, run_list=None, all_passes=True, ratdb_tag=None):
    """Generate the file list by production labels
    """
    # Now loop through each module and get the file lists
    if not module:
        module = get_module_list_label(data_type, label)
    if len(module)==0:
        print "No modules for label %s" % label
        sys.exit()
    names = []
    sizes = []
    guids = []
    adlers = []
    for i, mod in enumerate(module):
        print "Generating file lists, %d of %d" % (i, len(module))
        n, s, g, a = get_list_by_label_module(data_type, label, mod, run_range, run_list=run_list, all_passes=all_passes, ratdb_tag=ratdb_tag)
        names += n
        sizes += s
        guids += g
        adlers += a
    if len(names)==0:
        print "No files for label %s, modules %s" % (label, module)
        sys.exit()
    return names, sizes, guids, adlers


def generate_list_by_version(version, data_type, module, run_range, run_list=None, all_passes=True, ratdb_tag=None):
    """Generate the file list by rat version
    """
    # Now loop through each module and get the file lists
    if not module:
        module = get_module_list_version(data_type, version)
    if len(module)==0:
        print "No modules for rat %s" % version
        sys.exit()
    names = []
    sizes = []
    guids = []
    adlers = []
    for i, mod in enumerate(module):
        print "Generating file lists, %d of %d" % (i, len(module))
        n, s, g, a= get_list_by_module_version(data_type, mod, version, run_range, run_list=run_list, all_passes=all_passes, ratdb_tag=ratdb_tag)
        names += n
        sizes += s
        guids += g
        adlers += a
    if len(names)==0:
        print "No files for version %s, modules %s" % (version, module)
        sys.exit()
    return names, sizes, guids, adlers

def generate_list_by_module_run(data_type, module, run_range, run_list=None, all_passes=True, ratdb_tag=None):
    """Generate the file list by module and run
    """
    # Now loop through each module and get the file lists
    if not module:
        module = get_module_list_version(data_type)
    if len(module)==0:
        print "No modules for rat %s" % version
        sys.exit()
    names = []
    sizes = []
    guids = []
    adlers = []
    for i, mod in enumerate(module):
        print "Generating file lists, %d of %d" % (i, len(module))
        n, s, g, a= get_list_by_module(data_type, mod, run_range, run_list=run_list, all_passes=all_passes, ratdb_tag=ratdb_tag)
        names += n
        sizes += s
        guids += g
        adlers += a
    if len(names)==0:
        print "No files for modules %s" % module
        sys.exit()
    return names, sizes, guids, adlers
    
def generate_list(parser):
    '''Main function called by both production_list and processing_list
    '''
    global copy_type, server_list
    args = parser.parse_args()

    if args.file_type not in allowed_file_types:
        print "Unknown filetype, allowed types are:"
        print ", ".join(f for f in allowed_file_types)
        sys.exit()

    # First, check that no help commands are requested
    if args.list_labels:
        database.connect_db(args.db_server, args.db_port, args.db_name)
        for l in sorted(get_labels(args.file_type, args.module)):
            print l
        sys.exit()
    if args.list_versions:
        database.connect_db(args.db_server, args.db_port, args.db_name)
        for v in sorted(get_versions(args.file_type, args.module)):
            print v
        sys.exit()
    if args.list_modules:
        database.connect_db(args.db_server, args.db_port, args.db_name)
        for m in sorted(get_modules(args.file_type, args.label, args.version)):
            print m
        sys.exit()

    # Check that the output file is OK
    if os.path.exists(args.output):
        overwrite = raw_input("%s exists, overwrite [y/N]?: " % args.output)
        if overwrite!="y" and overwrite!="Y":
            print "Specify another filename"
            sys.exit()

    # Connect to the database
    database.connect_db(args.db_server, args.db_port, args.db_name)

    # If there is a run range, need exactly two arguments
    run_range = None
    run_list = None
    if args.__contains__('exactly') and args.exactly:
        try:
            with open(args.exactly) as f:
                exactly = [line.split('\t') for line in f]
                exactly = [(module,int(run),int(passno)) for module,run,passno in exactly]
        except Exception as e:
            print 'Could not parse the exactly file: %s' % str(e)
            sys.exit(1)
    elif args.run_range:
        run_range = args.run_range
        if len(args.run_range)!=2:
            print "Require two arguments for run range (lower upper)"
            sys.exit()
    elif args.__contains__('run_file') and args.run_file:
        with open(args.run_file[0]) as f:
            run_list = [int(line) for line in f]
        # The +1 is needed b/c couchdb query is non-inclusive over a range
        run_range = [min(run_list), max(run_list)+1]

    if grid.copy == grid.srm_copy or not args.__contains__('copy_type') or args.copy_type is None or args.copy_type is 'surl':
        # Use SRMs to copy from 
        copy_type = "SURL"
        server_list = grid.get_server_preferences()

    if args.__contains__('all_passes'):
        all_passes = args.all_passes
    else:
        all_passes = True
        
    if args.__contains__('ratdb_tag'):
        ratdb_tag = args.ratdb_tag
    else:
        ratdb_tag = None
        
    if args.__contains__('exactly') and args.exactly:
        if ratdb_tag:
            parser.error('Cannot use --exactly and --ratdb-tag together')
        names, sizes, guids, adlers = generate_list_from_exactly(args.file_type, exactly)
    elif args.label:
        names, sizes, guids, adlers = generate_list_by_label(args.label, args.file_type, args.module, run_range, run_list=run_list, all_passes=all_passes, ratdb_tag=ratdb_tag)
    elif args.version:
        names, sizes, guids, adlers = generate_list_by_version(args.version, args.file_type, args.module, run_range, run_list=run_list, all_passes=all_passes, ratdb_tag=ratdb_tag)
    else:
        names, sizes, guids, adlers = generate_list_by_module_run(args.file_type, args.module, run_range, run_list=run_list, all_passes=all_passes, ratdb_tag=ratdb_tag)

    # And write the output file
    utilities.write_grabber_file(args.output, copy_type, names, sizes, guids, adlers)
