package org.workflowsim.planning;

import java.util.*;
import org.cloudbus.cloudsim.Log;
import org.workflowsim.CondorVM;
import org.workflowsim.Task;
import java.time.*;


public class aco extends BasePlanningAlgorithm {

	public static int no_tasks;
	public static int no_vms;
	public static List<Task> cloudletList = new LinkedList<>();
	public static List<CondorVM> vmList = new LinkedList<>();
	public static HashMap <Task,HashMap<CondorVM,Double>> expected_exec_time = new HashMap <Task,HashMap<CondorVM,Double>>();
	public static int tmax = 100;
	public static double initialPheromone;
	public static double Q = 100;
	public static double alpha = 0.3;
	public static double beta = 1;
	public static double rho = 0.4;
	public static int ants = 4;
    public static double makespan = 0.0;
    public static double [] time = new double[100];
    
    public static List <CondorVM> schedule() {
		for(int i=0;i<no_tasks;i++){
			
			HashMap <CondorVM,Double> tmp = new HashMap <CondorVM,Double>(); 
			for (int j = 0; j < no_vms; j++) {
				double val = (double)cloudletList.get(i).getCloudletLength()/(vmList.get(j).getNumberOfPes()*vmList.get(j).getMips()) + (double)cloudletList.get(i).getCloudletLength()/vmList.get(j).getBw();
				tmp.put(vmList.get(j),val);
			}
			expected_exec_time.put(cloudletList.get(i),tmp);
		}
		List <CondorVM> allocatedtasks = new LinkedList<>();
		int stp = 0;
		for(int i=0;i<cloudletList.size();i+=no_vms){
			List <Task> tempCloudletList = new LinkedList<>();
			for(int j=i;j<cloudletList.size() && j<i+no_vms;j++){
				tempCloudletList.add(cloudletList.get(j));
			}
			List<CondorVM> temp_List = implementScheduling(tempCloudletList, vmList, tmax,stp);
			stp++;
			for(int j=0;j<temp_List.size();j++){
				allocatedtasks.add(temp_List.get(j));
			}
		}

		return allocatedtasks;
	}
    
    public static List <CondorVM> implementScheduling(List<Task> cloudletList, List<CondorVM> vmList, int tmax,int stp) {
		int cur_itr = 1;
		double opt_len=1000000.0;
		ArrayList <CondorVM> opt_sol = new ArrayList<>();
		HashMap <Task,HashMap<CondorVM,Double>> pheromone = new HashMap <Task,HashMap<CondorVM,Double>>();
		for (int i = 0; i < cloudletList.size(); i++) {
			HashMap <CondorVM,Double> tmp = new HashMap <CondorVM,Double>(); 
			for (int j = 0; j < no_vms; j++) {
				tmp.put(vmList.get(j),1.0);
			}
			pheromone.put(cloudletList.get(i),tmp);
		}
        
		while (cur_itr <= tmax) {
			Random rand = new Random();
			ArrayList<ArrayList<CondorVM>> tabu = new ArrayList<>();
			for (int i = 0; i < ants; i++) {
				tabu.add(new ArrayList<>());
				int randindex = rand.nextInt(no_vms);
				CondorVM randVM = vmList.get(randindex);
				tabu.get(i).add(randVM);
				for (int k = 1; k < cloudletList.size(); k++) {
					double mx = 0.0;
					int selected_vm = 0;
					for (int l = 0; l < no_vms; l++) {
						if (tabu.get(i).contains(vmList.get(l)) == false) {
							double pj = Math.pow(pheromone.get(cloudletList.get(k)).get(vmList.get(l)), alpha)
									* Math.pow(1.0 / (expected_exec_time.get(cloudletList.get(k)).get(vmList.get(l))),beta);
							if (pj > mx) {
								mx = pj;
								selected_vm = l;
							}
						}
					}
					tabu.get(i).add(vmList.get(selected_vm));
				}
			}

			double global_mn = 1000000.0;
			int selected_ant = 0;
			for (int i = 0; i < ants; i++) {
				double mx = 0.0;
				for (int j = 0; j < no_vms; j++) {
					double len = time[j];
					for (int k = 0; k < tabu.get(i).size(); k++) {
						if (tabu.get(i).get(k) == vmList.get(j)) {
							len += expected_exec_time.get(cloudletList.get(k)).get(vmList.get(j));
						}
					}
					if (mx < len) {
						mx = len;
					}
				}
				for (int j = 0; j < tabu.get(i).size(); j++) {
					CondorVM selected_vm = tabu.get(i).get(j);
					double delta_p = Q / mx;
					pheromone.get(cloudletList.get(j)).put(selected_vm,(1 - rho) * pheromone.get(cloudletList.get(j)).get(selected_vm) + delta_p);
				}
				if (global_mn > mx) {
					global_mn = mx;
					selected_ant = i;
				}
			}
			double delta_p = Q / global_mn;
			for (int i = 0; i < cloudletList.size(); i++) {
				CondorVM selected_vm = tabu.get(selected_ant).get(i);
				pheromone.get(cloudletList.get(i)).put(selected_vm, pheromone.get(cloudletList.get(i)).get(selected_vm) + delta_p);
			}
			if(global_mn < opt_len){
				opt_len = global_mn;
				opt_sol = tabu.get(selected_ant);
			}
			cur_itr++;
		}

        for(int j=0;j<no_vms;j++){
			double len=time[j];
			for(int k=0;k<cloudletList.size();k++){
				if(opt_sol.get(k)==vmList.get(j)){
					len+=expected_exec_time.get(cloudletList.get(k)).get(vmList.get(j));
				}
			}
			time[j] = len;
		}
		return opt_sol;
	}

   

    public void run() {
    	Instant start = Instant.now();
        Log.printLine("ACO running with " + getTaskList().size()
                + " tasks.");

        for (Object vmObject : getVmList()) {
            CondorVM vm = (CondorVM) vmObject;
            vmList.add(vm);
        }

          cloudletList = getTaskList();
          no_tasks = cloudletList.size();
          no_vms = vmList.size();
          List <CondorVM> scheduledTasks = schedule();
          
          Instant end = Instant.now();
          Duration timeElapsed = Duration.between(start, end);
          Log.printLine("\n\nExecution Time: "+ timeElapsed.getNano() +"\n");
          double makespan = 0;
          for(int i=0;i<no_vms;i++){
              makespan = Math.max(makespan, time[i]);
          }
 
          Log.printLine("\nMakespan: "+ makespan +"\n");
         for(int i=0;i<scheduledTasks.size();i++) {
       	  cloudletList.get(i).setVmId(scheduledTasks.get(i).getId());
       	  Log.printLine(cloudletList.get(i).getCloudletId() + " : " + cloudletList.get(i).getCloudletLength() + " : " + scheduledTasks.get(i).getId());
         }
         

        System.out.println("\n\n\n\n\n\n\n");
        
    }

    
}
