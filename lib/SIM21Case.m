classdef SIM21Case < handle
    properties
        name;
        dataPath;
        tmpDataPath;
        outputPath;
        ncube
        fstar;
        vbc;
        vc;
        fx;
        sed;
        tau;
        feedback;
        delayParam;
        pop;
        fsfunc;
        phVersion;
        zeta;
        ID;
    end

    methods
        function obj = SIM21Case(name,dataPath,tmpDataPath,outputPath,ncube,fstar,vbc,vc,fx,sed,tau,feedback,delayParam,pop,fsfunc,phVersion,zeta)
            obj.name = name;
            obj.dataPath = dataPath;
            obj.tmpDataPath = tmpDataPath;
            obj.outputPath = outputPath;
            obj.ncube = ncube;
            obj.fstar = fstar;
            obj.vbc = vbc;
            obj.vc = vc;
            obj.fx = fx;
            obj.sed = sed;
            obj.tau = tau;
            obj.feedback = feedback;
            obj.delayParam = delayParam;
            obj.pop = pop;
            obj.fsfunc = fsfunc;
            obj.phVersion = phVersion;
            obj.zeta = zeta;
            obj.ID = obj.getID();
        end


        function runID = getID(obj)
            runID = ['_',num2str(obj.ncube),'_',num2str(obj.fstar),'_',num2str(obj.vbc),'_',num2str(obj.vc),'_',num2str(obj.fx),...
                     '_',num2str(obj.sed),'_',num2str(obj.tau),'_',num2str(obj.feedback),'_',num2str(obj.delayParam),'_',num2str(obj.pop),...
                     '_',num2str(obj.fsfunc),'_',num2str(obj.phVersion)];
        end
    end
end