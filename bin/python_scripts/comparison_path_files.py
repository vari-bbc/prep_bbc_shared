import os
import hashlib
import logging

from collections import OrderedDict

logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d:%H:%M:%S',
                    level=logging.DEBUG)

# this is for turning off the logging
logging.getLogger().disabled = False


class ComparisonPathOfFiles(object):
    def __init__(self, BackupPathOfFiles, sourcePathOfFiles):
        # path of files in the previous directory
        self._backup_path_of_files = BackupPathOfFiles
        # path of files in the source directory
        self._source_path_of_files = sourcePathOfFiles
        # abs path used
        self._hardlinks_path_from_previous_dir_list = None
        # abs path used
        self._copy_files_path_from_source_dir_list = None

        # symlinks dict from backup_dir
        self._symlinks_dict_from_backup_dir = BackupPathOfFiles.symlink_dict
        # symlinks dict from source_dir
        self._symlinks_dict_from_source_dir = sourcePathOfFiles.symlink_dict

        self._symlinks_verified_from_source_dir = OrderedDict()
        self.compare_files()

    @property
    def symlinks_dict_from_backup_dir(self):
        return self._symlinks_dict_from_backup_dir

    @property
    def symlinks_dict_from_source_dir(self):
        return self._symlinks_dict_from_source_dir

    @property
    def backup_dir_prefix_path(self):
        return self._backup_path_of_files.prefix_path

    @property
    def source_dir_prefix_path(self):
        return self._source_path_of_files.prefix_path

    @property
    def backup_path_of_files(self):
        return self._backup_path_of_files

    @property
    def source_path_of_files(self):
        return self._source_path_of_files

    @property
    def hardlinks_path_from_previous_dir_list(self):
        return self._hardlinks_path_from_previous_dir_list

    @property
    def copy_files_path_from_source_dir_list(self):
        return self._copy_files_path_from_source_dir_list

    def __str__(self):
        return "hard-links for dest dir : {}\n copy files for dest dir {}".format(
            self.hardlinks_path_from_previous_dir_list,
            self.copy_files_path_from_source_dir_list)

    def compare_files(self):

        '''
        This function is for comparing only files. Not symlinks. setter function
        :return: None
        '''
        hardlinks_path_from_previous_dir_list = list()
        copy_files_path_from_source_dir_list = list()

        pre_relative_path_of_files_list = [i.replace(self.backup_dir_prefix_path, "") for i in
                                           self.backup_path_of_files.abs_path_of_files_list]

        for abs_path_src_file in self.source_path_of_files.abs_path_of_files_list:
            logging.info(abs_path_src_file)

            src_file = abs_path_src_file.replace(self.source_dir_prefix_path, "")
            if src_file in pre_relative_path_of_files_list:

                abs_pre_file = os.path.abspath(self.backup_dir_prefix_path + "/" + src_file)
                assert (os.path.isfile(abs_pre_file))

                abs_src_file = os.path.abspath(self.source_dir_prefix_path + "/" + src_file)
                assert (os.path.isfile(abs_pre_file))

                logging.info("---begin to calculate md5sum---")
                logging.info("This is an abs path of a file in the backup directory and its md5sum.")
                logging.debug(abs_pre_file)

                abs_pre_file_md5 = self.calculate_md5sum(abs_pre_file)
                logging.debug(abs_pre_file_md5)

                logging.info("This is an abs path of a file in the source directory and its md5sum.")
                logging.debug(abs_src_file)

                abs_src_file_md5 = self.calculate_md5sum(abs_src_file)
                logging.debug(abs_src_file_md5)
                logging.info('---End calculating md5sum---')
                logging.info("\n")

                if abs_src_file_md5 == abs_pre_file_md5:
                    hardlinks_path_from_previous_dir_list.append(abs_pre_file)
                ## if files are not identical
                else:
                    copy_files_path_from_source_dir_list.append(abs_src_file)
            ## if a file in source directory is not found in the previous directory
            else:
                ## abs_src_file = os.path.join(source_dir, src_file)
                copy_files_path_from_source_dir_list.append(abs_path_src_file)

        self._hardlinks_path_from_previous_dir_list = hardlinks_path_from_previous_dir_list
        self._copy_files_path_from_source_dir_list = copy_files_path_from_source_dir_list
        return None

    def calculate_md5sum(self, a_file):
        '''

        calculating md5sum
        :param a_file: an abs path of file
        :return: md5sum value
        '''
        with open(a_file, 'rb')as fin:
            file_hash = hashlib.md5()
            chunk = fin.read(8192)
            while chunk:
                file_hash.update(chunk)
                chunk = fin.read(8192)
        return file_hash.hexdigest()

    def compare_symlinks(self):
        '''
        This function is for comparing symlinks. Not hard links (unable to identify hard links in Python) - setter function
        :return: None
        '''

        ## calling for only a key which is abs_symlinks
        relative_path_backup_files = [i.replace(self.backup_dir_prefix_path) for i in
                                      self.symlinks_dict_from_backup_dir.keys()]
        for abs_source_symlinks_file in self.symlinks_dict_from_source_dir:

            ## creating relative symlinks string
            src_symlinks_file = abs_source_symlinks_file.replace(self.source_dir_prefix_path + "/", "")

            ## if relative src symlinks in symlinks in the backup directory
            ## Then, keep the symlinks in a dictionary
            if src_symlinks_file in relative_path_backup_files:
                self._symlinks_verified_from_source_dir[abs_source_symlinks_file] = self.symlinks_dict_from_source_dir[
                    abs_source_symlinks_file]

            else:
                logging.info("The below symlink file is not found in the backup dir.")
                logging.debug(abs_source_symlinks_file)
                self._symlinks_verified_from_source_dir[abs_source_symlinks_file] = self.symlinks_dict_from_source_dir[
                    abs_source_symlinks_file]