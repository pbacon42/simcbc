
from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import table as ligolw_table
from glue.ligolw import lsctables

class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
        pass

xmldoc = ligolw_utils.load_filename("test.xml",contenthandler=LIGOLWContentHandler, verbose = True)

l = xmldoc.getElements(lambda e: (len(e.tagName) > 1))

table = ligolw_table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
#coinc_table = ligolw_table.get_table(xmldoc, lsctables.CoincTable.tableName)
#coinc_def_table = ligolw_table.get_table(xmldoc, lsctables.CoincDefTable.tableName)
#coinc_map_table = ligolw_table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
