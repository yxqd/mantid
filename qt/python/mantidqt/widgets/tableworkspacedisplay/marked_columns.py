from mantidqt.widgets.tableworkspacedisplay.error_column import ErrorColumn


class MarkedColumns:
    X_LABEL = "[X{}]"
    Y_LABEL = "[Y{}]"
    Y_ERR_LABEL = "[Y{}_YErr]"

    def __init__(self):
        self.as_x = []
        self.as_y = []
        self.as_y_err = []  # type: list[ErrorColumn]
        # self.as_x_err = []

    def _add(self, col_index, add_to, remove_from):
        assert all(
            add_to is not remove for remove in remove_from), "Can't add and remove from the same list at the same time!"
        self._remove(col_index, remove_from)

        if col_index not in add_to:
            add_to.append(col_index)

    def _remove(self, col_index, remove_from):
        """
        Remove the column index from all lists
        :param col_index: The column index to be removed
        :type remove_from: list[list[Union[int, ErrorColumn]]]
        :param remove_from: List of lists from which the column index will be removed
        :return:
        """
        for list in remove_from:
            try:
                list.remove(col_index)
            except ValueError:
                # column not in this list, but might be in another one so we continue the loop
                continue

        # if the column previously had a Y Err associated with it -> this will remove it from the YErr list
        self._remove_associated_yerr_columns(col_index)

    def add_x(self, col_index):
        self._add(col_index, self.as_x, [self.as_y, self.as_y_err])

    def add_y(self, col_index):
        self._add(col_index, self.as_y, [self.as_x, self.as_y_err])

    def add_y_err(self, err_column):
        # remove all labels for the column index
        len_before_remove = len(self.as_y)
        self._remove(err_column, [self.as_x, self.as_y, self.as_y_err])

        # Check if the length of the list with columns marked Y has shrunk
        # -> This means that columns have been removed, and the label_index is now _wrong_
        # and has to be decremented to match the new label index correctly
        len_after_remove = len(self.as_y)
        if err_column.error_for_column > err_column.source_column and len_after_remove < len_before_remove:
            err_column.label_index -= (len_before_remove - len_after_remove)
        self.as_y_err.append(err_column)

    def remove_column(self, col_index):
        self._remove(col_index, [self.as_x, self.as_y, self.as_y_err])

    def _remove_associated_yerr_columns(self, col_index):
        # we can only have 1 Y Err for Y, so iterating and removing's iterator invalidation is not an
        # issue as the code will exit immediately after the removal
        for col in self.as_y_err:
            if col.error_for_column == col_index:
                self.as_y_err.remove(col)
                break

    def _make_labels(self, list, label):
        return [(col_num, label.format(index),) for index, col_num in enumerate(list)]

    def build_labels(self):
        extra_labels = []
        extra_labels.extend(self._make_labels(self.as_x, self.X_LABEL))
        extra_labels.extend(self._make_labels(self.as_y, self.Y_LABEL))
        err_labels = [(err_col.source_column, self.Y_ERR_LABEL.format(err_col.label_index),) for index, err_col in
                      enumerate(self.as_y_err)]
        extra_labels.extend(err_labels)
        return extra_labels

    def find_yerr(self, selected_columns):
        yerr_for_col = {}

        # for each selected column
        for col in selected_columns:
            # find the marked error column
            for yerr_col in self.as_y_err:
                # if found append the YErr's source column - so that the data from the columns
                # can be retrieved for plotting the errors
                if yerr_col.error_for_column == col:
                    yerr_for_col[col] = yerr_col.source_column

        return yerr_for_col
