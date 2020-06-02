AliAnaysisTaskHyperTriton3VtxPerf * AddMyTask ( TString name = " name ") {
  AliAnalysisManager * mgr = AliAnalysisManager::GetAnalysisManager();
  TString fileName = AliAnalysisManager :: GetCommonFileName();
  fileName += ": MyTask "; // create a subfolder in the file
  // now we create an instance of your task
  AliAnaysisTaskHyperTriton3VtxPerf * task = new AliAnaysisTaskHyperTriton3VtxPerf(true, name.Data());
  // add your task to the manager
  mgr->AddTask(task);
  // connect the manager to your task
  mgr->ConnectInput(task ,0 ,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task ,1 ,mgr->CreateContainer("MyOutputContainer ", TList::Class (), AliAnalysisManager::kOutputContainer, fileName.Data())) ;
  // important : return a pointer to your task
  return task ;
}