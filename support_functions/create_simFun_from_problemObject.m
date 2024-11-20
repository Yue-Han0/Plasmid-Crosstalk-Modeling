function F = create_simFun_from_problemObject(problemObject,additional_track_species)

    % Create a simFunction based on a problemObject 
    
    Mobj = problemObject.Model; 
    param_names = {problemObject.Estimated.Name};

    % Get observables from response map 
    response_map = char(problemObject.ResponseMap); 
    observable_raw = response_map(1:strfind(response_map,'=')-2); 
    observable_raw2 = strrep(observable_raw,'[',''); 
    observables = strrep(observable_raw2,']',''); 
    
    if exist("additional_track_species",'var')
        observables = [{observables} additional_track_species];
    end

    % Create dosed object for createSimFunction based on scheduled dose
    all_problemObject_doses = {problemObject.Doses.TargetName}; 
    all_problemObject_doses_reshape = reshape(all_problemObject_doses,size(problemObject.Doses)); 
    dosed = all_problemObject_doses_reshape(1,:); 

    % CreateSimFunction 
    F = createSimFunction(Mobj,param_names,observables,dosed); 

end