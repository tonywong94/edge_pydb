import h5py
from astropy.table import Table
from edge_pydb import util
from astropy import units as u

filename = util.fetch('NGC4047.2d_smo7.hdf5')
print('Working on',filename)
h5f = h5py.File(filename, 'r')
paths = [key for key in h5f.keys() if "__table_column_meta__" not in key]
print(paths)

# CO moments table
tab = Table.read(filename, path='comom_dil')
itab = tab.info(out=None)
itab['units'] = [u.Unit(key).to_string('latex_inline') for key in itab['unit']]
itab['colname']  = [key.replace('_','{\_}') for key in itab['name']]
itab.rename_column('n_bad', 'n{\_}bad')
itab['description'] = itab['description'].astype('U62')
# Rewrite some descriptions
replace_list = {"ra " : "R.A. ", 
                "dec " : "Decl. ",
                "radius" : "galactocentric radius",
                "azang" : "azimuthal angle",
                "13co " : "$^{13}$CO ",
                "co " : "CO ",
                "dil" : "dilated",
                "wgtd" : "weighted",
                "H2" : "H$_2$"}
for i_row, row in enumerate(itab['description']):
    for key in replace_list.keys():
        if key in row:
            row = row.replace(key,replace_list[key])
    itab['description'][i_row] = row        
# Finalize the table
itab = itab['colname','units','description','n{\_}bad']
itab.write('table_comom.tex' ,format='ascii.aastex',
           latexdict = {'tabletype': 'deluxetable*'},
           caption='Description of columns in the comom{\_}dil table.',
           overwrite=True)

# SSP table
tab = Table.read(filename, path='SSP_sm')
itab = tab.info(out=None)
for i_row, row in enumerate(itab['unit']):
    if row in ['none', 'fraction']:
        itab['unit'][i_row] = ''
itab['units'] = [u.Unit(key).to_string('latex_inline') for key in itab['unit']]
itab['colname']  = [key.replace('_','{\_}') for key in itab['name']]
itab['description'] = itab['description'].astype('U64')
itab.rename_column('n_bad', 'n{\_}bad')
itab['n{\_}bad'] = itab['n{\_}bad'].astype('U4')
# Rewrite some units
replace_list = {"1 \\times " : "", 
                "solMass / arcsec2" : "$\mathrm{M_{\odot}\,arcsec^{-2}}$"}
for i_row, row in enumerate(itab['units']):
    for key in replace_list.keys():
        if key in row:
            row = row.replace(key,replace_list[key])
    itab['units'][i_row] = row        
# Rewrite some descriptions
replace_list = {"ra " : "R.A. ", 
                "dec " : "Decl. ",
                "radius" : "galactocentric radius",
                "azang" : "azimuthal angle",
                "wgtd" : "weighted",
                "attnuation" : "attenuation",
                "H2" : "H$_2$",
                "<" : "$<$",
                "_" : "{\_}"}
for i_row, row in enumerate(itab['description']):
    for key in replace_list.keys():
        if key in row:
            row = row.replace(key,replace_list[key])
    itab['description'][i_row] = row        
# Finalize the table
itab = itab['colname','units','description','n{\_}bad']
itab.insert_row(3,['{\dots}','{\dots}','{\dots}','{\dots}'])
itab.write('table_ssp.tex' ,format='ascii.aastex',
           latexdict = {'tabletype': 'deluxetable*'},
           caption='Description of columns in the SSP table.',
           overwrite=True)


# flux_elines table
tab = Table.read(filename, path='flux_elines_sm')
itab = tab.info(out=None)
# Remove rows for lines except Halpha and Hbeta
exclude = ['[OIII]', '[OII]', '[OI]', '[NII]', '[SII]']
truth = [any(s in name for s in exclude) for name in itab['name']]
badrows = [i for i, x in enumerate(truth) if x]
itab.remove_rows(badrows)
# Remove common rows with RA, DEC, etc.
itab.remove_rows([3,4,5,6,7,8])
print(itab)
itab['units'] = [u.Unit(key).to_string('latex_inline') for key in itab['unit']]
itab['colname']  = [key.replace('_','{\_}') for key in itab['name']]
itab['description'] = itab['description'].astype('U64')
itab.rename_column('n_bad', 'n{\_}bad')
itab['n{\_}bad'] = itab['n{\_}bad'].astype('U4')
# Rewrite some units
replace_list = {"1 \\times " : "", 
                "solMass / arcsec2" : "$\mathrm{M_{\odot}\,arcsec^{-2}}$"}
for i_row, row in enumerate(itab['units']):
    for key in replace_list.keys():
        if key in row:
            row = row.replace(key,replace_list[key])
    itab['units'][i_row] = row        
# Rewrite some descriptions
replace_list = {"ra " : "R.A. ", 
                "dec " : "Decl. ",
                "radius" : "galactocentric radius",
                "azang" : "azimuthal angle",
                "wgtd" : "weighted",
                "attnuation" : "attenuation",
                "Halpha" : "H$\\alpha$",
                "Hbeta" : "H$\\beta$",
                "Ha " : "H$\\alpha$ ",
                "H2" : "H$_2$",
                "o3n2" : "O3N2",
                "<" : "$<$",
                ">" : "$>$",
                "_" : "{\_}"}
for i_row, row in enumerate(itab['description']):
    for key in replace_list.keys():
        if key in row:
            row = row.replace(key,replace_list[key])
    itab['description'][i_row] = row        
# Finalize the table
itab = itab['colname','units','description','n{\_}bad']
itab.insert_row(3,['{\dots}','{\dots}','{\dots}','{\dots}'])
itab.write('table_felines.tex' ,format='ascii.aastex',
           latexdict = {'tabletype': 'deluxetable*'},
           caption='Description of columns in the flux{\_}elines table.',
           overwrite=True)

