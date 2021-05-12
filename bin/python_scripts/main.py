import argparse
import logging

from comparison_path_files import ComparisonPathOfFiles
from copy_files_and_hard_links import CopyFilesHardlinks
from path_of_files import PathOfFiles

logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d:%H:%M:%S',
                    level=logging.DEBUG)

# turning off the logging function
logging.getLogger().disabled = False


def main(backup_dir, source_dir, dest_dir):
    backup_dir_files_path = PathOfFiles(backup_dir)
    source_dir_files_path = PathOfFiles(source_dir)
    comparison = ComparisonPathOfFiles(backup_dir_files_path, source_dir_files_path)
    start = CopyFilesHardlinks(comparison, dest_dir)
    start.run()


def test(pre_dir, source_dir, dest_dir):
    pre_dir_files_path = PathOfFiles(pre_dir)
    print(pre_dir_files_path.symlink_dict)

    logging.info("backup-dir-prefix_path")
    logging.info(pre_dir_files_path.prefix_path)

    source_dir_files_path = PathOfFiles(source_dir)
    logging.info("source-dir-prefix-path")
    logging.info(source_dir_files_path.prefix_path)

    tmp = ComparisonPathOfFiles(pre_dir_files_path, source_dir_files_path)

    print("hardlinks - unchanged")
    print(tmp.hardlinks_path_from_previous_dir_list)

    print("copy files - changed")
    print(tmp.copy_files_path_from_source_dir_list)
    start = CopyFilesHardlinks(tmp, dest_dir)
    start.run()

    logging.info("hard links backup to dest")
    logging.debug(start.hard_links_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Python scripts for doing the rsync link-dest job")

    parser.add_argument('-b',
                        '--backup_dir',
                        action='store',
                        type=str,
                        required=True,
                        help="backup-dir",
                        )

    parser.add_argument('-s',
                        '--source_dir',
                        action='store',
                        type=str,
                        required=True,
                        help='source-dir',
                        )
    parser.add_argument('-d',
                        '--dest_dir',
                        action='store',
                        type=str,
                        required=True,
                        help='dest-dir',
                        )

    args = parser.parse_args()

    main(args.backup_dir, args.source_dir, args.dest_dir)
