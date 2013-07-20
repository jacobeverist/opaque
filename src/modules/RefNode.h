

class RefNode {
public:

	double refX;
	double refZ;
	double refP;
	double gndX;
	double gndZ;
	double gndP;

	int nodeID;
	int jointID;
	int numJoints;

	double nom_ang[100];

	double maxError;
	int active_time;
	double errSum;
	double currError;
	int maxErrorReached;
	int newRoot;


    RefNode(int nid, int jid, double x, double z, double p, int numJoints, double *joints);
    ~RefNode();

	double getRefPoseX();
	double getRefPoseZ();
	double getRefPoseP();

	int getNodeID();
	int getJointID();
	void setGroundTruthPose(double x, double z, double p);

	double getGroundTruthPoseX();
	double getGroundTruthPoseZ();
	double getGroundTruthPoseP();

	double computeStabilityError(double *joints);
	double getStabilityError();
	double getMaxStabilityError();
	double getAvgStabilityError();
	int isMaxErrorReached();
	void updateTime();
	void setNewRoot(int isRoot);
	int isNewRoot();
};


