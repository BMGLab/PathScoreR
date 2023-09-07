
#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <queue>

using namespace Rcpp;
using namespace std;
vector<int> dfs_max(vector<vector<int>>* edges);
vector<vector<int>> bfs_ndv(vector<vector<int>>* edges,int max_d_len);
int dfs_len(vector<vector<int>>* edges, unordered_set<int>* visited,unsigned int edge);
vector<double> ndv_mean_calc(vector<vector<int>>* ndv, vector<int>* max_d,int max_d_len);
vector<vector<double>> P_numerator(vector<double>* ndv_mean,vector<vector<double>>* ndv_real,int max_d_len);
vector<vector<double>> ndv_real_val(vector<vector<int>>* ndv,vector<int>* max_d,int max_d_len);
vector<vector<double>> P_delimetor(vector<double>* ndv_mean,vector<vector<double>>* ndv_real,int max_d_len);
vector<vector<double>> P_ij(vector<vector<double>>* P_num,vector<vector<double>>* P_del);
vector<vector<int>> bfs_d(vector<vector<int>>* edges);
vector<double> S_i_back(const CharacterVector& rowNames,const CharacterVector& edges_1,
                        const CharacterVector& edges_2,vector<double>* p_values,
                        const double alpha, const int D);
vector<double> L_i_back(const vector<double>* p_values);
double log_sum_back(const vector<double>* p_values);
vector<double> sum_L_S(vector<double>* L_values,vector<double>* S_values);
vector<vector<double>> calculate_means(vector<vector<vector<double>>>* group_items,
                                       int gene_count,int group_count, int one_group_len);
vector<vector<double>> calculate_SD(vector<vector<vector<double>>>* group_items,
                                    int gene_count,int group_count, int one_group_len);
vector<vector<double>> calculate_SE(vector<vector<vector<double>>>* group_items,
                                    int gene_count,int group_count, int one_group_len,int sample_count);
DataFrame calculate_anova(vector<vector<vector<double>>>* group_items, // anova hesaplanacak // gene_count , group_count, one_group_len
                          int gene_count,int group_count, int one_group_len,int sample_count);
DataFrame calculate_t_stat(vector<vector<vector<double>>>* group_items, // t_stat hesaplanacak // gene_count , group_count, one_group_len
                           int gene_count,int group_count, int one_group_len,int sample_count);
vector<string> col_names(int gene_count);
// [[Rcpp::export]]
DoubleVector log_sum(const List& v) {
    vector<double> vec;
    vec.assign( v.begin(),v.end());
    unsigned int len = vec.size();
    double raw_total = 0.0;
    for(unsigned int i = 0;i<len;i++){
      raw_total += -log(vec.at(i));
    }
    DoubleVector total = DoubleVector::create(raw_total);
    return total;
}
double log_sum_back(const vector<double>* p_values) {
  unsigned int len = p_values->size();
  double raw_total = 0.0;
  for(unsigned int i = 0;i<len;i++){
    raw_total += -log(p_values->at(i));
  }
  return raw_total;
}
// [[Rcpp::export]]
DoubleVector L_i(const List& v){
  vector<double> vec;
  vec.assign(v.begin(),v.end());
  DoubleVector total = log_sum(v);
  double denominator = total.at(0);
  unsigned int len = v.size();
  DoubleVector result;
  for(unsigned int i = 0;i<len;i++){
    result.push_back(-log(vec.at(i))/denominator);
  }
  return result;
}
vector<double> L_i_back(const vector<double>* p_values){
  double total = log_sum_back(p_values);
  unsigned int len = p_values->size();
  vector<double> result;
  for(unsigned int i = 0;i<len;i++){
    result.push_back(-log(p_values->at(i))/total);
  }
  return result;
}
DataFrame calculate_p_values(vector<vector<vector<double>>>* group_items,
                                  int gene_count,int group_count, int one_group_len,int sample_count){
  if(group_count == 2){
    return calculate_t_stat(group_items,gene_count,group_count,one_group_len,sample_count);
  }
  return calculate_anova(group_items,gene_count,group_count,one_group_len,sample_count);
  
}
DataFrame calculate_anova(vector<vector<vector<double>>>* group_items, // gene_count , group_count, one_group_len
                                int gene_count,int group_count, int one_group_len,int sample_count){
  DataFrame df;
  for(unsigned int i = 0;i<gene_count;i++){
    vector<vector<double>> tmp;
    tmp.assign(group_items->at(i).begin(),group_items->at(i).end());
    DoubleVector tmp_values;
    CharacterVector tmp_groups;
    
    for(unsigned int j = 0;j<group_count;j++){
      for(unsigned int k = 0;k<one_group_len;k++){
        tmp_values.push_back(tmp.at(j).at(k));
        tmp_groups.push_back(char(48+j));
      }
    }
    df.push_back(tmp_values);
    df.push_back(tmp_groups);
  }
  return df;
  
}
vector<string> col_names(int gene_count){
  vector<string> cols;
  for(int i = 0;i<gene_count;i++){
    string tmp_1 = "";
    tmp_1 += char(48 + i);
    string tmp_2 = "";
    tmp_2 = char(48 + i) + 'g';
    cols.push_back(tmp_1);
    cols.push_back(tmp_2);
  }
  return cols;
}
DataFrame calculate_t_stat(vector<vector<vector<double>>>* group_items, // t_stat hesaplanacak // gene_count , group_count, one_group_len
                                int gene_count,int group_count, int one_group_len,int sample_count){
  DataFrame df;
  for(unsigned int i = 0;i<gene_count;i++){
    vector<vector<double>> tmp;
    tmp.assign(group_items->at(i).begin(),group_items->at(i).end());
    
    for(unsigned int j = 0;j<group_count;j++){
      DoubleVector tmp_values;
      for(unsigned int k = 0;k<one_group_len;k++){
        tmp_values.push_back(tmp.at(j).at(k));
      }
      df.push_back(tmp_values);
    }
    
  }
  return df;
}
vector<vector<double>> calculate_SE(vector<vector<vector<double>>>* group_items,
                            int gene_count,int group_count, int one_group_len,int sample_count){
  vector<vector<double>>&& sd_values = calculate_SD(group_items,gene_count,group_count,one_group_len);
  vector<vector<double>> result;
  for(unsigned int i = 0;i<gene_count;i++){
    vector<double> tmp;
    for(unsigned int j = 0;j<group_count;j++){
      tmp.push_back(sd_values.at(i).at(j)/sqrt(sample_count));
    }
  }
  return result;
}

vector<vector<double>> calculate_SD(vector<vector<vector<double>>>* group_items,
                            int gene_count,int group_count, int one_group_len){
  vector<vector<double>>&& group_means = calculate_means(group_items,gene_count,group_count,one_group_len);
  vector<vector<double>> result;
  for(unsigned int i = 0;i<gene_count;i++){
    vector<double> tmp;
    for(unsigned int j = 0;j<group_count;j++){
      double total = 0.0;
      for(unsigned k = 0;k<one_group_len;k++){
        total += pow(group_items->at(i)[j][k]-group_means.at(i).at(j),2);
      }
      tmp.push_back(sqrt(total/(one_group_len-1)));
    }
    result.push_back(tmp);
  }
  return result;
}
vector<vector<double>> calculate_means(vector<vector<vector<double>>>* group_items,
                                       int gene_count,int group_count, int one_group_len){
  vector<vector<double>> result;
  for(unsigned int i = 0;i<gene_count;i++){
    vector<double> tmp;
    for(unsigned int j = 0;j<group_count;j++){
      double total = 0.0;
      for(unsigned k = 0;k<one_group_len;k++){
        total += group_items->at(i)[j][k];
      }
      tmp.push_back(total/one_group_len);
    }
    result.push_back(tmp);
  }
  return result;
}
// [[Rcpp::export]]
DataFrame run_2(const List& list,const DoubleVector& p,const CharacterVector& rowNames,
                const CharacterVector& edges_1,const CharacterVector& edges_2,
                const NumericVector& groups, const double alpha, const int D){
  vector<double> p_values;
  p_values.assign(p.begin(),p.end());
  vector<int> group_vec;
  group_vec.assign(groups.begin(),groups.end());
  int group_count = 0;
  for(unsigned int i = 0;i<group_vec.size();i++){
    group_count = max(group_count,group_vec.at(i));
    group_vec[i] -= 1;
  }
  //cout<<"1"<<endl;
  vector<double>&& S_i_values = S_i_back(rowNames,edges_1,edges_2,&p_values,alpha,D);
  //cout<<"2"<<endl;
  vector<double>&& L_i_values = L_i_back(&p_values);
  //cout<<"3"<<endl;
  vector<double>&& sum_LS = sum_L_S(&L_i_values,&S_i_values);
  //cout<<"4"<<endl;
  vector<vector<double>> T_groups;
  for(unsigned int i = 0;i<group_count;i++){
    vector<double> T_i;
    for(unsigned int j = 0;j<list.size();j++){
      if(i == group_vec.at(j)){
        DoubleVector tmp = list.at(j);
        double total = 0.0;
        for(unsigned int k = 0;k<tmp.size();k++){
          total += tmp.at(k)*sum_LS.at(k);
        }
        T_i.push_back(total);
      }
    }
    T_groups.push_back(T_i);
  }
  //cout<<"5"<<endl;
  DataFrame result;
  for(unsigned int i = 0;i<T_groups.size();i++){
    DoubleVector tmp;
    for(unsigned int j = 0;j<list.size()/group_count;j++){
      tmp.push_back(T_groups.at(i).at(j));
    }
    result.push_back(tmp);
  }
  //cout<<"6"<<endl;
  return result;
}
// [[Rcpp::export]]
DataFrame run_1(const List& list,NumericVector& groups,const CharacterVector& rowNames,
                 const CharacterVector& edges_1,const CharacterVector& edges_2){
  vector<int> group_vec;
  int group_count = 0;
  group_vec.assign(groups.begin(),groups.end());
  for(unsigned int i = 0;i<group_vec.size();i++){
    group_count = max(group_count,group_vec.at(i));
    group_vec[i]-=1;
  }
  int one_group_len = group_vec.size()/group_count;
  vector<vector<vector<double>>> group_items; // gene_count , group_count, one_group_len
  for(unsigned int i = 0;i<rowNames.size();i++){
    vector<vector<double>> tmp_i;
    for(unsigned int j = 0;j<group_count;j++){
      vector<double> tmp_j;
      for(unsigned int k = 0;k<list.size();k++){
        DoubleVector tmp_k = list.at(k);
        if(group_vec.at(k) == j){
          tmp_j.push_back(tmp_k.at(i));
        }
      }
      tmp_i.push_back(tmp_j);
    }
    group_items.push_back(tmp_i);
  }
  return calculate_p_values(&group_items,rowNames.size(),group_count,one_group_len,list.size());
  
}
vector<double> sum_L_S(vector<double>* L_values,vector<double>* S_values){
  vector<double> result;
  for(unsigned int i = 0;i<L_values->size();i++){
    result.push_back(L_values->at(i)+S_values->at(i));
  }
  return result;
}
vector<double> S_i_back(const CharacterVector& rowNames,const CharacterVector& edges_1,
                 const CharacterVector& edges_2,vector<double>* p_values,
                 const double alpha, const int D){
  unordered_map<string,unsigned int> indexes;
  vector<string> row_names;
  vector<string> edge_1_names;
  vector<string> edge_2_names;
  edge_2_names.assign(edges_2.begin(),edges_2.end());
  edge_1_names.assign(edges_1.begin(),edges_1.end());
  row_names.assign(rowNames.begin(),rowNames.end());
  //cout<<"1.1"<<endl;
  unsigned int len = row_names.size();
  //cout<< len << endl;
  for(unsigned int i = 0;i<len;i++){
    //cout<< row_names.at(i)<<endl;
    indexes[row_names.at(i)] = i;
  }
  vector<vector<int>> edge_matrix;
  //cout<<"1.2"<<endl;
  for(unsigned int i = 0;i<len;i++){
    vector<int> tmp;
    for(unsigned int j = 0;j<len; j++){
      tmp.push_back(0);
    }
    edge_matrix.push_back(tmp);
  }
  //cout<<"1.3"<<endl;
  unsigned int len_edge_1 = edges_1.size();
  for(unsigned int i = 0;i<len_edge_1;i++){
    int col = indexes[edge_1_names.at(i)];
    int row = indexes[edge_2_names.at(i)];
    edge_matrix[col][row] = 1;
    edge_matrix[row][col] = 1;
  }/*
 for(unsigned int i = 0;i<edge_matrix.size();i++){
 for(unsigned int j = 0;j<edge_matrix.size();j++){
 cout<<edge_matrix[i][j]<< " ";
 }
 cout<<endl;
 }*/
  vector<int> max_d = dfs_max(&edge_matrix);
  //cout<<"1.4"<<endl;
  /*
   for(unsigned int i = 0;i<max_d.size();i++){
   cout<<max_d[i]<< endl;
   }*/
  vector<vector<int>>&& d_matrix = bfs_d(&edge_matrix); 
  //cout<<"1.5"<<endl;
  int max_d_len = 0;
  for(unsigned int i = 0;i<edge_matrix.size();i++){
    for(unsigned int j = 0;j<edge_matrix.size();j++){
      max_d_len = max(max_d_len,d_matrix.at(i).at(j));
    }
  }
  //cout<<"1.6"<<endl;
  vector<vector<int>>&& ndv = bfs_ndv(&edge_matrix,max_d_len);
  /*
   for(unsigned int i = 0;i<edge_matrix.size();i++){
   for(unsigned int j = 0;j<max_d_len;j++){
   cout<<ndv[i][j]<< " ";
   }
   cout<<endl;
   }
   */
  //cout<<max_d_len<<endl;
  
  
  vector<double>&& ndv_mean = ndv_mean_calc(&ndv,&max_d,max_d_len);
  /*
   for(unsigned int i = 0;i<ndv_mean.size();i++){
   cout<< ndv_mean[i]<<endl;
   }
   */
  vector<vector<double>>&& ndv_real = ndv_real_val(&ndv,&max_d,max_d_len);
  /*
   for(unsigned int i = 0;i<ndv_real.size();i++){
   for(unsigned int j = 0;j< max_d_len;j++){
   cout<<ndv_real[i][j]<< " "; 
   }
   cout<<endl;
   }
   */
  vector<vector<double>>&& P_num = P_numerator(&ndv_mean,&ndv_real,max_d_len);
  /*
   for(unsigned int i = 0;i<ndv.size();i++){
   for(unsigned int j = 0;j< ndv.size();j++){
   cout<<P_num[i][j]<< " "; 
   }
   cout<<endl;
   }
   */
  vector<vector<double>>&& P_del = P_delimetor(&ndv_mean,&ndv_real,max_d_len);
  /*
   for(unsigned int i = 0;i<ndv.size();i++){
   for(unsigned int j = 0;j< ndv.size();j++){
   cout<<P_del[i][j]<< " "; 
   }
   cout<<endl;
   }
   */
  vector<vector<double>>&& P_ij_matrix = P_ij(&P_num,&P_del);
  /*
   for(unsigned int i = 0;i<ndv.size();i++){
   for(unsigned int j = 0;j< ndv.size();j++){
   cout<<P_ij_matrix[i][j]<< " "; 
   }
   cout<<endl;
   }
   */
  //cout<<"1.7"<<endl;
  //cout<<P_ij_matrix.size()<<endl;
  //cout<<d_matrix.size()<<endl;
  //cout<<p_values->size()<<endl;
  double total = 0.0;
  for(unsigned int i = 0;i<d_matrix.size();i++){
    for(unsigned int j = 0;j<d_matrix.size();j++){
      if(d_matrix.at(i).at(j) != 0 && d_matrix.at(i).at(j)<D && p_values->at(i)<alpha){
        total += (fabs(P_ij_matrix.at(i).at(j))/d_matrix.at(i).at(j));
      }
    }
  }
  //cout<<"1.8"<<endl;
  vector<double> result;
  for(unsigned int i = 0;i<d_matrix.size();i++){
    double tmp_total = 0.0;
    for(unsigned int j = 0;j<d_matrix.size();j++){
      if(d_matrix.at(i).at(j)!= 0 && d_matrix.at(i).at(j)<D && p_values->at(i)<alpha){
        tmp_total += (fabs(P_ij_matrix.at(i).at(j))/d_matrix.at(i).at(j));
      }
    }
    result.push_back(tmp_total/total);
  }
  return result;
}

// [[Rcpp::export]]
DoubleVector S_i(const CharacterVector& rowNames,const CharacterVector& edges_1,
                 const CharacterVector& edges_2, const double alpha, const int D){
  unordered_map<string,unsigned int> indexes;
  vector<string> row_names;
  vector<string> edge_1_names;
  vector<string> edge_2_names;
  edge_2_names.assign(edges_2.begin(),edges_2.end());
  edge_1_names.assign(edges_1.begin(),edges_1.end());
  row_names.assign(rowNames.begin(),rowNames.end());
  unsigned int len = row_names.size();
  //cout<< len << endl;
  for(unsigned int i = 0;i<len;i++){
    //cout<< row_names.at(i)<<endl;
    indexes[row_names.at(i)] = i;
  }
  vector<vector<int>> edge_matrix;
  for(unsigned int i = 0;i<len;i++){
    vector<int> tmp;
    for(unsigned int j = 0;j<len; j++){
      tmp.push_back(0);
    }
    edge_matrix.push_back(tmp);
  }
  unsigned int len_edge_1 = edges_1.size();
  for(unsigned int i = 0;i<len_edge_1;i++){
    int col = indexes[edge_1_names.at(i)];
    int row = indexes[edge_2_names.at(i)];
    //cout<<col<<" "<<row<<endl;
    edge_matrix[col][row] = 1;
    edge_matrix[row][col] = 1;
  }/*
  for(unsigned int i = 0;i<edge_matrix.size();i++){
    for(unsigned int j = 0;j<edge_matrix.size();j++){
      cout<<edge_matrix[i][j]<< " ";
    }
    cout<<endl;
  }*/
  vector<int>&& max_d = dfs_max(&edge_matrix);
  /*
  for(unsigned int i = 0;i<max_d.size();i++){
    cout<<max_d[i]<< endl;
  }*/
  vector<vector<int>>&& d_matrix = bfs_d(&edge_matrix); 
  int max_d_len = 0;
  for(unsigned int i = 0;i<edge_matrix.size();i++){
    for(unsigned int j = 0;j<edge_matrix.size();j++){
      //cout<<d_matrix.at(i).at(j)<<" ";
      max_d_len = max(max_d_len,d_matrix.at(i).at(j));
    }
    //cout<<endl;
  }
  //cout<<max_d_len<<endl;
  vector<vector<int>>&& ndv = bfs_ndv(&edge_matrix,max_d_len);
  for(unsigned int i = 0;i<edge_matrix.size();i++){
    for(unsigned int j = 0;j<max_d_len;j++){
      cout<<ndv[i][j]<< " ";
    }
    cout<<endl;
  }
  //cout<<max_d_len<<endl;
  
  
  vector<double>&& ndv_mean = ndv_mean_calc(&ndv,&max_d,max_d_len);
  /*
  for(unsigned int i = 0;i<ndv_mean.size();i++){
    cout<< ndv_mean[i]<<endl;
  }
   */
  vector<vector<double>>&& ndv_real = ndv_real_val(&ndv,&max_d,max_d_len);
  /*
  for(unsigned int i = 0;i<ndv_real.size();i++){
    for(unsigned int j = 0;j< max_d_len;j++){
      cout<<ndv_real[i][j]<< " "; 
    }
    cout<<endl;
  }
   */
  vector<vector<double>>&& P_num = P_numerator(&ndv_mean,&ndv_real,max_d_len);
  /*
  for(unsigned int i = 0;i<ndv.size();i++){
    for(unsigned int j = 0;j< ndv.size();j++){
      cout<<P_num[i][j]<< " "; 
    }
    cout<<endl;
  }
   */
  vector<vector<double>>&& P_del = P_delimetor(&ndv_mean,&ndv_real,max_d_len);
  /*
  for(unsigned int i = 0;i<ndv.size();i++){
    for(unsigned int j = 0;j< ndv.size();j++){
      cout<<P_del[i][j]<< " "; 
    }
    cout<<endl;
  }
   */
  vector<vector<double>>&& P_ij_matrix = P_ij(&P_num,&P_del);
  /*
  for(unsigned int i = 0;i<ndv.size();i++){
    for(unsigned int j = 0;j< ndv.size();j++){
      cout<<P_ij_matrix[i][j]<< " "; 
    }
    cout<<endl;
  }
   */
  
  double total = 0.0;
  for(unsigned int i = 0;i<d_matrix.size();i++){
    for(unsigned int j = 0;j<d_matrix.size();j++){
      if(d_matrix.at(i).at(j) != 0 && d_matrix.at(i).at(j)<D){
        total += (fabs(P_ij_matrix.at(i).at(j))/d_matrix.at(i).at(j));
      }
    }
  }
  DoubleVector result;
  for(unsigned int i = 0;i<d_matrix.size();i++){
    double tmp_total = 0.0;
    for(unsigned int j = 0;j<d_matrix.size();j++){
      if(d_matrix.at(i).at(j)!= 0 && d_matrix.at(i).at(j)<D){
        tmp_total += (fabs(P_ij_matrix.at(i).at(j))/d_matrix.at(i).at(j));
      }
    }
    result.push_back(tmp_total/total);
  }
  //cout<<ndv.size()<<endl;
  /*
 
  for(unsigned int i = 0;i<ndv.size();i++){
    for(unsigned int j = 0;j< max_d_len;j++){
      cout<<ndv[i][j]<< " "; 
    }
    cout<<endl;
  }*/
  //cout<<"done"<<endl;
  /*
  for(unsigned int i = 0;i<ndv.size();i++){
    NumericVector tmp;
    for(unsigned int j = 0;j < max_d_len;j++){
      tmp.push_back(ndv.at(i).at(j));
    }
    result.push_back(tmp);
  }
   */
  return result;
}
vector<vector<double>> P_ij(vector<vector<double>>* P_num,vector<vector<double>>* P_del){
  vector<vector<double>> result;
  for(unsigned int i = 0;i<P_num->size();i++){
    vector<double> tmp; 
    for(unsigned int j = 0;j<P_del->size();j++){
      if(i==j){
        tmp.push_back(1);
        continue;
      }
      tmp.push_back(P_num->at(i)[j]/P_del->at(i)[j]);
    }
    result.push_back(tmp);
  }
  return result;
}

vector<vector<double>> P_delimetor(vector<double>* ndv_mean,vector<vector<double>>* ndv_real,int max_d_len){
  vector<vector<double>> result;
  for(unsigned int i = 0;i<ndv_mean->size();i++){
    vector<double> tmp;
    
    for(unsigned int j = 0;j<ndv_mean->size();j++){
      if(i==j){
        tmp.push_back(1);
        continue;
      }
      double total_1 = 0.0;
      double total_2 = 0.0;
      for(unsigned int k = 0;k<max_d_len;k++){
        total_1 += pow(ndv_real->at(j)[k]-ndv_mean->at(j),2);
        total_2 += pow(ndv_real->at(i)[k]-ndv_mean->at(i),2);
      }
      double total = sqrt(total_1)*sqrt(total_2);
      if(total == 0){
        total = 1;
      }
      tmp.push_back(total);
    }
    result.push_back(tmp);
  }
  return result;
}
vector<vector<int>> bfs_d(vector<vector<int>>* edges){
  queue<int> q;
  vector<vector<int>> result;
  for(unsigned int i = 0;i<edges->size();i++){
    vector<int> tmp;
    for(unsigned int j = 0;j<edges->size();j++){
      tmp.push_back(0);
    }
    result.push_back(tmp);
  }
  for(int i = 0;i<edges->size();i++){
    q.push(i);
    unordered_set<int> visited;
    vector<int> one_line;
    int level = 0;
    while(!q.empty()){
      int current = q.back();
      //cout << current << endl;
      level++;
      q.pop();
      vector<int> tmp;
      for(int j = 0;j< edges->size();j++){
        if(edges->at(current).at(j) == 1 && visited.find(j) == visited.end()){
          q.push(j);
          visited.insert(j);
          tmp.push_back(j);
        }
      }
      for(unsigned int j = 0;j<tmp.size();j++){
        //cout<<level<<endl;
        result.at(i).at(tmp.at(j)) = level;
      }
    }
  }
  return result;
}

vector<vector<double>> P_numerator(vector<double>* ndv_mean,vector<vector<double>>* ndv_real,int max_d_len){
  vector<vector<double>> result;
  for(unsigned int i = 0;i<ndv_mean->size();i++){
    vector<double> tmp;
    
    for(unsigned int j = 0;j<ndv_mean->size();j++){
      if(i==j){
        tmp.push_back(1);
        continue;
      }
      double total = 0.0;
      for(unsigned int k = 0;k<max_d_len;k++){
        total += (double)(ndv_real->at(i)[k]-ndv_mean->at(i))*(ndv_real->at(j)[k]-ndv_mean->at(j));
      }
      tmp.push_back(total);
    }
    result.push_back(tmp);
  }
  return result;
}
vector<vector<double>> ndv_real_val(vector<vector<int>>* ndv,vector<int>* max_d,int max_d_len){
  vector<vector<double>> ndv_real;
  for(unsigned int i = 0;i<ndv->size();i++){
    vector<double> tmp;
    for(unsigned int j = 0;j<max_d_len;j++){
      tmp.push_back(((double)ndv->at(i)[j])/((double)max_d->at(i)));
    }
    ndv_real.push_back(tmp);
  }
  return ndv_real;
}
vector<double> ndv_mean_calc(vector<vector<int>>* ndv, vector<int>* max_d,int max_d_len){
  vector<double> ndv_mean;
  for(unsigned int i = 0;i<ndv->size();i++){
    double total = 0.0;
    for(unsigned int j = 0;j<max_d_len;j++){
      total += ndv->at(i)[j];
    }
    total /= (double) max_d->at(i);
    ndv_mean.push_back(total);
  }
  return ndv_mean;
}
vector<int> dfs_max(vector<vector<int>>* edges){
  vector<int> tmp;
  for(unsigned int i = 0;i<edges->size();i++){
    unordered_set<int> visited;
    tmp.push_back(dfs_len(edges,&visited,i));
  }
  return tmp;
}
int dfs_len(vector<vector<int>>* edges, unordered_set<int>* visited,unsigned int edge){
  int tmp = 0;
  for(unsigned int i = 0;i<edges->size();i++){
    if(edges->at(edge).at(i)==0||i == edge || visited->find(edge) != visited->end()){
      continue;
    }
    visited->insert(i);
    tmp = max(tmp,dfs_len(edges,visited,i));
    visited->erase(visited->find(i));
  }
  return 1 + tmp;
}
vector<vector<int>> bfs_ndv(vector<vector<int>>* edges,int max_d_len){
  queue<int> q;
  vector<vector<int>> result;
  for(int i = 0;i<edges->size();i++){
    q.push(i);
    unordered_set<int> visited;
    vector<int> one_line;
    while(!q.empty()){
      int current = q.back();
      //cout << current << endl;
      q.pop();
      int level = 0;
      for(int j = 0;j< edges->size();j++){
        if(edges->at(current).at(j) == 1 && visited.find(j) == visited.end()){
          level++;
          //cout << j << endl;
          q.push(j);
          visited.insert(j);
        }
      }
      //cout<< level<< endl;
      one_line.push_back(level);
    }
    while(one_line.size() < max_d_len){
      one_line.push_back(0);
    }
    result.push_back(one_line);
  }
  return result;
}
/*** R
col_names_gene_count <- function(gene_count){
  cols <- c()
  for(i in 1:gene_count){
    cols[2*i-1] = as.character(i)
    cols[2*i] = paste(as.character(i),"g",sep = "")
  }
  return(cols)
}
calculate_p_values_t_stat <- function(gene_count,datas){
  p_values <- c()
  for(i in 1:gene_count){
    tmp_names_1<- paste(as.character(i),"g",sep = "")
    tmp_p_val = t.test(datas[[as.character(i)]] ,datas[[tmp_names_1]], var.equal=TRUE)$p.value
    p_values[i] = as.double(tmp_p_val)
  }
  return(p_values)
}
calculate_p_values_anova <- function(gene_count,datas){
  p_values <- c()
  for(i in 1:gene_count){
    tmp_names_1<- paste(as.character(i),"g",sep = "")
    tmp_names_1<- paste("`",tmp_names_1,sep = "")
    tmp_names_1<- paste(tmp_names_1,"`",sep = "")
    tmp_names_2<- paste("`",as.character(i),sep = "")
    tmp_names_2<- paste(tmp_names_2,"`",sep = "")
    tmp_names_2<- paste(tmp_names_2," ~ ",sep = "")
    tmp_names <- paste(tmp_names_2,tmp_names_1,sep = "")
    tmp_formula = as.formula(tmp_names)
    tmp_p_val<-summary(aov(formula = tmp_formula ,data = datas))
    tmp_p_val2<-as.list(tmp_p_val[[1]])
    p_values[i]=as.double(tmp_p_val2$`Pr(>F)`[1])
  }
  return(p_values)
}
col_names_groups <- function(group_count){
  cols <- c()
  for(i in 1:group_count){
    cols[i]=as.character(i)
  }
  return(cols)
}
run <- function(list,groups,row_names,edges_1,edges_2,alpha,D){
  tmp<-run_1(list,groups,row_names,edges_1,edges_2)
  group_count<- max(groups)
  cols <- col_names_gene_count(length(row_names))
  colnames(tmp) <- cols
  if(group_count == 2){
    t_test <- calculate_p_values_t_stat(length(row_names),tmp)
    final <- run_2(list,t_test,row_names,edges_1,edges_2,groups,alpha,D)
    colnames(final)<- col_names_groups(2)
    return(final)
  }else{
    anova <- calculate_p_values_anova(length(row_names),tmp)
    final <- run_2(list,anova,row_names,edges_1,edges_2,groups,alpha,D)
    colnames(final)<- col_names_groups(group_count)
    
    return(final)
  }
}
*/
