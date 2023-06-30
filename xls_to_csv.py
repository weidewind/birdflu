import xlrd
import csv
import re
import optparse

parser = optparse.OptionParser()
parser.add_option('--xls', help='', type='str')
parser.add_option('--csv', help='', type='str')

options, args = parser.parse_args()


wb = xlrd.open_workbook(options.xls)
# print(wb.sheet_names())
# sh = wb.sheet_by_index(0)
your_csv_file = open(options.csv, 'w')
wr = csv.writer(your_csv_file, delimiter='	', quoting=csv.QUOTE_MINIMAL)
start = 0
for sheet in wb.sheets():
	for rownum in range(start, sheet.nrows):
		vals = sheet.row_values(rownum)
		newvals = [re.sub('\n', ' ', str(s)) for s in vals]
		wr.writerow(newvals)
	start = 1

your_csv_file.close()