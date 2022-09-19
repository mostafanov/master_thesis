classdef employeea
    properties
        Name
        Category = 'trainer'
        IDnumber 
    end
    methods
        function E= set.Name (E,Name)
            %E single row car
            if ischar(Name) ,, ndims(Name)==2 ,, size (Name,1)==1 %ischar(Name)is character  ndims(Name)==2 two dimensional,, size (Name,1)==1 raw vector
                E.Name= Name; 
            else
                error ('invaled name');
            end 

        end
                function E= set.Category (E,Newcategory)
                    possCategory= ('Manager','trainer','CEO');
                    switch Newcategory
           
                        case possCategory
                E.Category= Newcategory;
                        otherwise
            
                error ('invaled Category');
            end 

                end
                function E= set.IDnumber (E,IDN)
            %E single row car
            if isnumeric(IDN) ,, isscalar(IDN) ,, cell(IDN)==floor(IDN) %floor greater than zero
                E.IDnumber= IDN; 
            else
                error ('invaled IDN');
            end 

        end
    end
end
