from openpyxl.utils import column_index_from_string

def delete_column(ws, delete_column):
    """
    This method empties a specified column from the worksheet. However,
    it still preserves the number of columns (ie. it doesn't actually
    delete any columns, but it empties one of the columns)

    :param ws:                  (Object) Worksheet (not Workbook)
    :param delete_column:       (String) Column Letter OR (Integer) Column Number
    :return:                    (Object) Worksheet
    """
    if isinstance(delete_column, str):
        delete_column = column_index_from_string(delete_column)
    assert delete_column >= 1, "Column numbers must be 1 or greater"

    for column in range(delete_column, ws.max_column + 1):
        for row in range(1, ws.max_row + 1):
            ws.cell(row=row, column=column).value = ws.cell(row=row, column=column+1).value

    return ws