import terminaltables as tt


class Tableist():

    def __init__(self, rows=1, cols=1, blank='---', style='single'):
        super().__init__()
        self.blank = blank
        header = [blank for i in range(cols)]
        self._grid = tt.AsciiTable([])
        [self._grid.table_data.append(header.copy()) for _ in range(rows)]
        self.style = style

    def __getitem__(self, row):
        if isinstance(row, slice):  # If slice, convert to list
            start = row.start or 0
            stop = row.stop or self.nrows
            step = row.step or 1
            row = [r for r in range(start, stop, step)]
        if isinstance(row, list):  # If list, process each location
            temp = [self._grid.table_data[r] for r in row]
            temp2 = Tableist(len(temp), self.ncols, style=self.style)
            temp2._grid.table_data = temp
            temp2._grid.inner_heading_row_border = False
            return temp
        return self._grid.table_data[row]

    # def __setitem__(self, row, cols):
    #     if isinstance(row, int):
    #         self._grid.table_data[row] = cols
    #     if isinstance(row, slice):  # If slice, convert to list
    #         start = row.start or 0
    #         stop = row.stop or self.nrows
    #         step = row.step or 1
    #         row = [r for r in range(start, stop, step)]
    #     if isinstance(row, list):


    def set_col(self, col, vals):
        try:
            len(vals)
        except TypeError:
            vals = [vals for r in range(self.nrows)]
        for r in range(self.nrows):
            self._grid.table_data[r][col] = vals[r]

    def set_row(self, row, vals):
        try:
            len(vals)
        except TypeError:
            vals = [vals for c in range(self.ncols)]
        for c in range(self.ncols):
            self._grid.table_data[row][c] = vals[c]

    def _set_style(self, style):
        if hasattr(self, '_grid'):
            data = self._grid.table_data
        else:
            data = []
        if style.lower().startswith('a'):
            self._style = 'ascii'
            self._grid = tt.AsciiTable(data)
        elif style.lower().startswith('d'):
            self._style = 'double'
            self._grid = tt.DoubleTable(data)
        elif style.lower().startswith('s'):
            self._style = 'single'
            self._grid = tt.SingleTable(data)
        elif style.lower()[0] in ['g', 'm']:
            self._style = 'markdown'
            self._grid = tt.GithubFlavoredMarkdownTable(data)
        self._grid.inner_heading_row_border = True
        self._grid.padding_left = 3
        self._grid.padding_right = 3
        self._grid.justify_columns = {col: 'center' for col in range(self.ncols)}

    def _get_style(self):
        return self._style

    style = property(fset=_set_style, fget=_get_style)

    def _set_blank(self, blank):
        if hasattr(self, '_blank'):
            for i, row in enumerate(self._grid.table_data):
                for j, col in enumerate(self._grid.table_data[i]):
                    if self._grid.table_data[i][j] == self._blank:
                        self._grid.table_data[i][j] = blank
        self._blank = blank

    def _get_blank(self):
        if hasattr(self, '_blank'):
            return self._blank

    blank = property(fset=_set_blank, fget=_get_blank)

    @property
    def ncols(self):
        header = [self.blank for i in self._grid.table_data[0]]
        return len(header)

    @property
    def nrows(self):
        return len(self._grid.table_data)

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def size(self):
        return self.nrows * self.ncols

    def get_row(self, row):
        if not isinstance(row, int):
            row = self.get_col(0)._grid.table_data.index([row])
        temp = self._grid.table_data[row]
        temp2 = Tableist(1, self.ncols)
        for col, item in enumerate(temp):
            temp2[0][col] = item
        return temp2

    def get_col(self, col):
        if not isinstance(col, int):
            col = self._grid.table_data[0].index(col)
        temp = [self._grid.table_data[row][col] for row in range(self.nrows)]
        temp2 = Tableist(self.nrows, 1)
        for row, item in enumerate(temp):
            temp2[row][0] = item
        return temp2

    def add_row(self, row=None, num=1):
        header = [self.blank for i in self._grid.table_data[0]]
        if row:
            self._grid.table_data.insert(row, header)
        else:
            self._grid.table_data.append(header)
        if num > 1:
            self.add_row(row, num=num-1)

    def add_col(self, col=None, num=1):
        for row in self._grid.table_data:
            if col:
                row.insert(col, self.blank)
            else:
                row.append(self.blank)
        if num > 1:
            self.add_col(col=col, num=num-1)

    def drop_col(self, col):
        for row in self._grid.table_data:
            row.pop(col)

    def drop_row(self, row):
        self._grid.table_data.pop(row)

    @property
    def index(self):
        return self.get_col(0)

    @property
    def header(self):
        return self.get_row(0)

    def __repr__(self):
        return self._grid.table.__str__()
