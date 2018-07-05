classdef superHandleClass < matlab.mixin.SetGet & dynamicprops
%     This is the abstract interface class for all Classes
%     it contains routines that apply to all objects that inherit from this
%     class, like: contructor; copy; clone; save; load; struc2obj;
%     obj2struc; delete
%     The HGSETGET class is an abstract class that provides an HG-style
%     property set and get interface.  HGSETGET is a subclass of HANDLE, so 
%     any classes derived from HGSETGET are handle classes.
%     (set,get,setdisp,getdisp)
%     The DYNAMICPROPS class is an abstract class that provides support
%     for dynamic properties for MATLAB objects.  Dynamic properties can
%     be used to attach temporary data to MATLAB objects.
%     (addprop)
    properties (Dependent = true)
        Type;% This Property enables to detect an object class when this latter has been converted to a structure;
    end
    methods %Constructor
        function obj = superHandleClass(varargin)
            if nargin==0,return;end% making an empty object
            switch class(varargin{1})
                case class(obj)
                    copy(varargin{1},obj);
                case 'struct'
                    struc2obj(obj,varargin{1});
                otherwise
                    % no structure or object given
            end
            if nargin==1, return; end% only one structure or object where given
            mobj = metaclass(obj);
            for p=1:1:length(mobj.Properties)% Check any Given Property/Value pair
                Input = find(strcmp(varargin, mobj.Properties{p}.Name));if~isempty(Input), obj.(mobj.Properties{p}.Name) = varargin{Input+1}; end
            end
        end
    end
    methods% Special Setting & Getting
        function out = get.Type(obj)
            out = class(obj);
        end
    end
    methods% Interface Functions
        function delete(obj) % delete object and all the objects in its properties
            mobj = metaclass(obj);
            % No delete required for: Abstract, Dependent,private and Transient Properties
            sel = find(cellfun(@(cProp)(~cProp.Constant && ... 
                ~cProp.Abstract && ...
                ~cProp.Transient && ...
                ~cProp.Dependent && ...
                ~strcmp(cProp.GetAccess,'private') && ...
                ~strcmp(cProp.GetAccess,'private')),mobj.Properties));
            if isempty(sel), return; end
            for p=sel(:)'
                try
                if ~iscell(obj.(mobj.Properties{p}.Name))
                   if obj.isH(obj.(mobj.Properties{p}.Name))&&~strcmp('Parent',mobj.Properties{p}.Name)% Dealing with non Parent Handle property
                      if isvalid(obj.(mobj.Properties{p}.Name)), delete(obj.(mobj.Properties{p}.Name)); end% Needs to be deleted      
                   end
                else
                    for l=1:1:length(obj.(mobj.Properties{p}.Name))
                        if obj.isH(obj.(mobj.Properties{p}.Name){l})&&~strcmp('Parent',mobj.Properties{p}.Name)% Dealing with non Parent Handle property
                           if isvalid(obj.(mobj.Properties{p}.Name){l}), delete(obj.(mobj.Properties{p}.Name){l});end% Needs to be deleted    
                        end
                    end
                end
                catch
                end
            end
            clear obj;
        end
        function cobj = clone(obj) % make an independent clone of the object
            cobj = eval(class(obj));
            copy(obj,cobj);
        end
        function copy(obj,cobj) % Transfer properties from one object to another
            if~strcmp(obj.Type,cobj.Type), error('Copy requires objects of the same class'); end
            mobj = metaclass(obj);
            % Do not copy: Abstract, Dependent,private and Transient Properties
            sel = find(cellfun(@(cProp)(~cProp.Constant && ... 
                ~cProp.Abstract && ...
                ~cProp.Transient && ...
                ~cProp.Dependent && ...
                ~strcmp(cProp.GetAccess,'private') && ...
                ~strcmp(cProp.GetAccess,'private')),mobj.Properties));
            for p=sel(:)'
                %disp([class(obj) ' : ' mobj.Properties{p}.Name]);
                if ~iscell(obj.(mobj.Properties{p}.Name))
                   if obj.isH(obj.(mobj.Properties{p}.Name))&&~strcmp('Parent',mobj.Properties{p}.Name)% Dealing with non Parent Handle property
                      cobj.(mobj.Properties{p}.Name) = clone(obj.(mobj.Properties{p}.Name));% Needs to be cloned      
                   else
                      cobj.(mobj.Properties{p}.Name) = obj.(mobj.Properties{p}.Name);
                   end
                else
                   prop = cell(1,length(obj.(mobj.Properties{p}.Name)));
                   for l=1:1:length(obj.(mobj.Properties{p}.Name))
                       if obj.isH(obj.(mobj.Properties{p}.Name){l})&&~strcmp('Parent',mobj.Properties{p}.Name)% Dealing with non Parent Handle property
                          prop{l} = clone(obj.(mobj.Properties{p}.Name){l});% Needs to be cloned      
                       else
                          prop{l} = obj.(mobj.Properties{p}.Name){l};
                       end
                   end
                   cobj.(mobj.Properties{p}.Name) = prop;
                end                
            end
        end
        function struc = saveobj(obj)
            struc = obj2struc(obj);
        end
        function struc = obj2struc(obj)%needed to save objects
            mobj = metaclass(obj);
            %disp('obj2struc');
            struc.Type = obj.Type;
            % Do not convert: Abstract, Dependent,private and Transient Properties
            sel = find(cellfun(@(cProp)(~cProp.Constant && ... 
                ~cProp.Abstract && ...
                ~cProp.Transient && ...
                ~cProp.Dependent && ...
                ~strcmp(cProp.GetAccess,'private') && ...
                ~strcmp(cProp.GetAccess,'private')),mobj.Properties));
            for p=sel(:)'
                %disp([class(obj) ' : ' mobj.Properties{p}.Name]);
                if ~iscell(obj.(mobj.Properties{p}.Name))
                   if obj.isH(obj.(mobj.Properties{p}.Name))&&~strcmp('Parent',mobj.Properties{p}.Name)% Dealing with non Parent Handle property
                      struc.(mobj.Properties{p}.Name) = obj2struc(obj.(mobj.Properties{p}.Name));% Needs to be converted      
                   elseif obj.isH(obj.(mobj.Properties{p}.Name))&&strcmp('Parent',mobj.Properties{p}.Name)% Dealing with Parent Handle Property
                      % do nothing
                   else
                      struc.(mobj.Properties{p}.Name) = obj.(mobj.Properties{p}.Name);
                   end
                else
                   prop = cell(1,length(obj.(mobj.Properties{p}.Name)));
                   for l=1:1:length(obj.(mobj.Properties{p}.Name))
                       if obj.isH(obj.(mobj.Properties{p}.Name){l})&&~strcmp('Parent',mobj.Properties{p}.Name)% Dealing with non Parent Handle property
                          prop{l} = obj2struc(obj.(mobj.Properties{p}.Name){l});% Needs to be converted      
                       elseif obj.isH(obj.(mobj.Properties{p}.Name){l})&&strcmp('Parent',mobj.Properties{p}.Name)% Dealing with Parent Handle Property
                          % do nothing
                       else
                          prop{l} = obj.(mobj.Properties{p}.Name){l};
                       end
                   end
                   struc.(mobj.Properties{p}.Name) = prop;
                end                
            end
        end
        function obj = struc2obj(obj,struc)%needed to load objects
            %disp('struc2obj');
            Properties = fieldnames(struc);
            for p=1:1:length(Properties)
                try % to put every property from the struct into a correct property of the object
                    if ~iscell(struc.(Properties{p}))
                        if isfield(struc.(Properties{p}),'Type')% need to convert to object
                           obj.(Properties{p}) = struc2obj(eval(struc.(Properties{p}).Type),struc.(Properties{p}));
                        else
                           obj.(Properties{p}) = struc.(Properties{p});
                        end
                    else
                        prop = cell(1,length(struc.(Properties{p})));
                        for l=1:1:length(struc.(Properties{p}))
                            if isfield(struc.(Properties{p}){l},'Type')% need to convert to object
                                prop{l} = struc2obj(eval(struc.(Properties{p}){l}.Type),struc.(Properties{p}){l});
                            else
                                prop{l} = struc.(Properties{p}){l};
                            end
                        end
                        obj.(Properties{p}) = prop;
                    end
                catch %#ok<CTCH>
                    %disp([Properties{p} ' failed to be loaded']);
                end
            end
        end
        function out = validH(obj,handle)
            % Property Field to test (handle) must be given as a String
            H = obj.(handle);% store in variable to eliminate risk for recursions
            out = superHandleClass.isH(H);% isH is a function defined in objectclass Utilities            
        end
    end
    methods (Static = true)
        function obj = loadobj(struc)
            obj = struc2obj(eval(struc.Type),struc);
        end
        function out = isH(in)
            % test whether in is a valid handle
            if isempty(in), out = false; return; end% empty check
            if ~isa(in,'handle'), out = false; return; end% handle check
            if ~isvalid(in), out = false; return; end% check whether handle has been deleted 
            out = true;% succeeded all tests
        end
        function out = convertUInt(in) % convert to the most efficient number of bits
           if isempty(in), out = []; return; end
           M = max(in(:));
           if M<=intmax('uint8'),out = uint8(in);return;end
           if M<=intmax('uint16'), out = uint16(in);return;end
           if M<=intmax('uint32'), out = uint32(in);return;end
           out = uint64(in);
       end
    end
end % classdef