import os
from collections import OrderedDict


class PathOfFiles(object):
    def __init__(self, pre_fix_path):
        self._prefix_path = os.path.abspath(pre_fix_path) + "/"
        self._abs_path_of_files_list = None
        self._symlink_dict = OrderedDict()
        self.get_path_of_files()

    @property
    def prefix_path(self):
        return self._prefix_path

    @property
    def symlink_dict(self):
        return self._symlink_dict

    @property
    def abs_path_of_files_list(self):
        return self._abs_path_of_files_list

    def get_path_of_files(self):
        '''
        get the absolute path of files from the previous directory or the source directory
        :set: set symlink_dict, abs_path_of_files_list
        '''
        os.chdir(self.prefix_path)
        abs_path_of_files_list = list()
        for abs_dir, sub_dirs, files in os.walk(self.prefix_path):

            if ".git" in abs_dir in abs_dir: #or ".snakemake" in abs_dir:
                pass
            else:
                for f in files:
                    tmp_path = os.path.join(abs_dir, f)
                    abs_file_path = os.path.abspath(tmp_path)

                    if os.path.islink(abs_file_path):
                        symlink_path_of_file = os.path.abspath(abs_file_path)
                        symlink_origin_path_of_file = os.readlink(abs_file_path)
                        self._symlink_dict[symlink_path_of_file] = symlink_origin_path_of_file

                    elif os.path.isfile(abs_file_path):
                        abs_path_of_files_list.append(abs_file_path)

                    else:
                        print(self.prefix_path, f)
                        raise ValueError("This is not a file nor a link")
        self._abs_path_of_files_list = abs_path_of_files_list

    def validate_files_path(self):
        '''
        validate files in lists
        :return: None or Error raised.
        '''

        if not os.path.isdir(self.prefix_path):
            print("This is a directory!")
            print(self.prefix_path)
            raise ValueError()
