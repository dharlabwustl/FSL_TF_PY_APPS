#!/bin/bash
######################
######## USER DETAILS#######
project_ID=${1}





# Function to get the column number given the column name
get_column_number() {
    local csv_file="$1"   # The CSV file to search in
    local column_name="$2" # The column name to find

    # Get the header (first line) of the CSV file
    header=$(head -n 1 "$csv_file")

    # Split the header into an array of column names
    IFS=',' read -r -a columns <<< "$header"

    # Loop through the columns and find the index of the column name
    for i in "${!columns[@]}"; do
        if [[ "${columns[$i]}" == "$column_name" ]]; then
            # Print the 1-based index (cut and other tools expect 1-based indexes)
            echo $((i  ))
            return
        fi
    done

    # If the column is not found, print an error and return a failure status
    echo "Column '$column_name' not found!" >&2
    return 1
}
curl  -u   $XNAT_USER:$XNAT_PASS  -X GET   $XNAT_HOST/data/projects/${project_ID}/experiments/?format=csv  > sessions.csv
sessions_list='sessions.csv' ##'sessions_COLI_ANALYTICS_STEP3_20231122041129_ordered.csv' ##'wrong_data_from_arjun.csv' ##vns_list_to_fix_df_1.csv' ###'sessions_COLI_ANALYTICS_STEP3_20231122041129_ordered.csv' 
csv_file=${sessions_list} 
column_name="ID" 
column_number=$(get_column_number "$csv_file" "$column_name")
echo ${column_number}
while IFS=',' read -ra array; do
    SESSION_ID=${array[${column_number}]}
    echo ${SESSION_ID}  
#    directory_to_create_destroy
#    ./STEPS_COMPARTMENTS_CSF.sh  ${SESSION_ID} ##${}SNIPR02_E09040 ##  SNIPR01_E00011  ##SNIPR02_E02991 ## ${SESSION_ID}
    break
done < <(tail -n +2 "${sessions_list}")
