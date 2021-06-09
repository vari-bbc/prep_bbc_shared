import logging
import os.path
import shutil

logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d:%H:%M:%S',
                    level=logging.DEBUG)

logging.getLogger().disabled = False


class CopyFilesHardlinks(object):
    def __init__(self, ComparisonPathFiles, path_of_dest_dir):
        self._backup_dir_prefix_path = ComparisonPathFiles.backup_dir_prefix_path
        self._source_dir_prefix_path = ComparisonPathFiles.source_dir_prefix_path

        # abs path of files in the previous directory in list
        self._hard_links_list = ComparisonPathFiles.hardlinks_path_from_previous_dir_list
        self._copy_files_path_from_source_dir_list = ComparisonPathFiles.copy_files_path_from_source_dir_list
        self._prefix_dest_dir_path = path_of_dest_dir + "/"
        self._symlinks_dict_from_source_dir = ComparisonPathFiles.symlinks_dict_from_source_dir

    def run(self):
        self.copy_changed_files()
        self.hard_links_unchanged_files()
        self.create_symlinks()

    @property
    def prefix_dest_dir_path(self):
        return self._prefix_dest_dir_path

    @property
    def hard_links_list(self):
        return self._hard_links_list

    @property
    def backup_dir_prefix_path(self):
        return self._backup_dir_prefix_path

    @property
    def source_dir_prefix_path(self):
        return self._source_dir_prefix_path

    def hard_links_unchanged_files(self):
        '''
        create hard links to unchanged files in the destination directory.
        :return: None
        '''
        # hard_links_list = abs_path_of_files_from_previous_dir
        for backup_dir_files in self._hard_links_list:

            # This is for getting something like this.
            # /secondary/projects/bbc/research/foobar_index/ecoli/sequence/e.coli.fa
            # ecoli/sequence/e.coli.fa
            relative_path_dest_file = backup_dir_files.replace(self.backup_dir_prefix_path, "")
            abs_path_dest_file = self.prefix_dest_dir_path + relative_path_dest_file

            # To create dirs
            full_path_dest_file = os.path.dirname(abs_path_dest_file)
            os.makedirs(full_path_dest_file, exist_ok=True)

            logging.info("---creating hard links to unchanged files---")
            logging.info(backup_dir_files + " -> " + abs_path_dest_file)
            # To create hard links to unchanged file
            os.link(backup_dir_files, abs_path_dest_file)

        logging.info("---creating hard links is done.---")

    def copy_changed_files(self):
        '''
        copy changed files from the source directory into the destination directory.
        :return: None
        '''

        for source_dir_files in self._copy_files_path_from_source_dir_list:
            relative_path_dest_file = source_dir_files.replace(self.source_dir_prefix_path, "")
            abs_path_dest_file = self.prefix_dest_dir_path + relative_path_dest_file

            # To create dirs
            full_path_dest_file = os.path.dirname(abs_path_dest_file)
            os.makedirs(full_path_dest_file, exist_ok=True)

            shutil.copy2(source_dir_files, abs_path_dest_file)

    def create_symlinks(self):
        '''
        create symlinks in the destination directory
        :return: None
        '''
        for abs_src_symlinks in self._symlinks_dict_from_source_dir:
            # This is a relative symlinks (file) path
            src_symlinks = abs_src_symlinks.replace(self.source_dir_prefix_path, "")

            change_dir = self.prefix_dest_dir_path + "/" + os.path.dirname(src_symlinks)
            logging.info(change_dir)

            os.makedirs(change_dir, exist_ok=True)
            os.chdir(change_dir)

            # Key is a relative path of the original file
            os.symlink(self._symlinks_dict_from_source_dir[abs_src_symlinks],
                       os.path.basename(src_symlinks))

            logging.info("---creating symlinks")
            logging.info('origin of symlinks file')
            logging.debug(self._symlinks_dict_from_source_dir[abs_src_symlinks])

            logging.info('symlinks')
            logging.debug(os.path.basename(src_symlinks))
