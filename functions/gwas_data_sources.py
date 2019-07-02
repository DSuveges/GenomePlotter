import json
import os.path
import sys

def get_db_properties(db_properties = "DEV2", db_list = False):
    '''
    Return database connection properties for spotrel (release snapshot database).
    '''

    # Path to json file:
    # databasePropertyFile = '/nfs/spot/sw/prod/gwas/scripts/ftpSummaryStatsRelease/databases.json'
    databasePropertyFile = '/Users/dsuveges/Project/GenomePlotter/databases.json'

    # test if file exists:
    if not os.path.isfile(databasePropertyFile):
        print("[Error] gwas database property file does not exist. Exiting.")
        # sys.exit()

    # Reading json with the credentials
    try:
        with open(databasePropertyFile) as json_file:  
            databaseProperties = json.load(json_file)
    except IOError:
        print("[Error] Failed to read gwas database property file. Exiting.")

    # Checking if only a list is expected:
    if db_list:
        db_instances = databaseProperties.keys().join(", ")
        print("[Info] The available database instances: " + db_instances)


    # Curation app database
    spotpro_properties_file = '/nfs/spot/sw/prod/gwas/config/application.properties'

    # Data release database
    spotrel_properties_file = '/nfs/spot/sw/prod/gwas/config/application_datarelease.properties'

    # DEV3 database properties - for testing script
    dev3_properties_file = '/nfs/spot/sw/prod/gwas/scripts/ftpSummaryStatsRelease/dev3_database.properties'

    if db_properties == 'SPOTPRO':
        properties_file = spotpro_properties_file
    elif db_properties == 'SPOTREL':
        properties_file = spotrel_properties_file
    elif db_properties == 'DEV3':
        properties_file = dev3_properties_file
    else:
        pass

    try:
        file = open(properties_file, 'r')
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
    else:
        for line in  file:
            if 'spring.datasource.url' in line:
                dsn_tns_property = line.split()
                conn, dsn_tns = dsn_tns_property[1].split('@')
                ip, port, sid = dsn_tns.split(':')

            if 'spring.datasource.username' in line:
                username_property = line.split()
                username = username_property[1]

            if 'spring.datasource.password' in line:
                password_property = line.split()
        file.close()

        return ip, port, sid, username, password


if __name__ == '__main__':
    db_config = 'dev3'
    get_db_properties(db_config)
