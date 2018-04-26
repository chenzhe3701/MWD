% chenzhe, 2018-03-18
% try superpixels() function, using exx data

addChenFunction;
load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt\WE43_T6_C1_s_all_grain_1144_local_map.mat', 'data');
exx=data(5).exx_local;
exx(isnan(exx))=0;

N = round(length(exx(:))/100);
[L,NumLabels] = superpixels(exx,N);
gb=find_boundary_from_ID_matrix(L);
myplot(exx,gb);


disp('required number of superpixels: ');
disp(N);
disp('actual number of superpixels: ');
disp(NumLabels);
