"""
    util.py

    Author: H.C. v. Stockhausen < hc at vst.io >
    Date:   2012-10-14
    
"""

def calc_table_col_widths(rows, table_width):
    max_chars_per_col = [0] * len(rows[0])
    for row in rows:
        for idx, col in enumerate(row):
            for line in str(col).split('\n'):
                max_chars_per_col[idx] = max(len(line),
                    max_chars_per_col[idx])
    sum_chars = sum(max_chars_per_col)
    return [(x * table_width / sum_chars) for x in max_chars_per_col]
 
