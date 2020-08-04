

#	./bin/2d_acoustic_rtm_pml_cpu_optimal_stack.e  ./par/parameter_2d_acoustic_rtm_optimal_stack_beijing\

#time mpirun	 -np  45	-machinefile nodefile  ./bin/2d_acoustic_rtm_pml_cpu_optimal_stack.e 	./par/parameter_2d_acoustic_rtm_waxian_v1   </dev/null>./log/log_opt  2>./log/err_opt &

#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_datalsm_layer_v1  </dev/null>log_layone  2>err_layone &


#time mpirun	 -np	4	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_haiyou_modify  </dev/null>log_a  2>err_layone &
#time mpirun	 -np	4	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_haiyou_modify_after_smooth  </dev/null>log_a  2>err_layone &
#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_no_water_ve1  </dev/null>log_a  2>err_layone &
#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_refraction_v3  </dev/null>log_b  2>err_layone &
#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_refracted_near_offset_v1  </dev/null>log_a  2>err_layone &
#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_refracted_all_offset_v1  </dev/null>log_b  2>err_layone &
#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_refracted_near_offset_v2  </dev/null>log_c  2>err_layone &
#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_refracted_all_offset_v2  </dev/null>log_d  2>err_layone &
#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_2801_641_diving_wave  </dev/null>log_a  2>err_layone &
#time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_no_water_diving_ve1  </dev/null>log_a  2>err_layone &
time mpirun	 -np	2	-machinefile nodefile   ./bin/2d_acoustic_rtm_pml_cpu_v1.e  ./par/parameter_2d_acoustic_rtm_no_water_diving_add_reflection_ve1  </dev/null>log_a  2>err_layone &


