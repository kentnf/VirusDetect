
#perl virus_detect.pl --thread_num 24 --coverage_cutoff 0.3 --depth_cutoff 10 KLL97
#perl virus_detect.pl --host_reference databases/tomato --thread_num 24 --coverage_cutoff 0.3 --depth_cutoff 10 MX > report_MX.txt &

perl virus_detect.pl --host_reference databases/tomato --thread_num 24 --coverage_cutoff 0.3 --depth_cutoff 10 --min_identity 97 test_data  > repor_test1.txt
perl virus_detect.pl --host_reference databases/tomato --thread_num 24 --coverage_cutoff 0.3 --depth_cutoff 10 --min_identity 80 test_data2 > repor_test2.txt

