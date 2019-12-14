//class definitons
class jets
{
    double pt, eta, phi, m;
    int btag;

    public:
    jets() = default;
    jets(double, double, double, double, int);

    void setAll(double, double, double, double, int);
    void setPt(double);
    void setEta(double);
    void setPhi(double);
    void setM(double);
    void setBTag(int);


    double getPt();
    double getEta();
    double getPhi();
    double getM();
    int    getBTag();
    TLorentzVector get4P();
};

// constructor and method definitions
inline jets::jets(double Pt, double Eta, double Phi, double M, int BTag)
{
    pt   = Pt;
    eta  = Eta;
    phi  = Phi;
    m    = M;
    btag = BTag;
};

void jets::setAll(double Pt, double Eta, double Phi, double M, int BTag)
{
    pt   = Pt;
    eta  = Eta;
    phi  = Phi;
    m    = M;
    btag = BTag;
};

void jets::setPt(double Pt){pt=Pt;}
void jets::setEta(double Eta){eta=Eta;}
void jets::setPhi(double Phi){phi=Phi;}
void jets::setM(double M){m=M;}
void jets::setBTag(int BTag){btag=BTag;}

double jets::getPt(){return pt;}
double jets::getEta(){return eta;}
double jets::getPhi(){return phi;}
double jets::getM(){return m;}
int    jets::getBTag(){return  btag;}

TLorentzVector jets::get4P()
{
    TLorentzVector dummy4P;
    dummy4P.SetPtEtaPhiM(pt, eta, phi, m);
    return dummy4P;
}

//class definitons
class leptons
{
    double pt, eta, phi, m, q;
    int pdgid;
    
    public:
    leptons() = default;
    leptons(double, double, double, double, double, int);

    void setAll(double, double, double, double, double, int);
    void setPt(double);
    void setEta(double);
    void setPhi(double);
    void setM(double);
    void setQ(double);
    void setPDGID(int);

    double getPt();
    double getEta();
    double getPhi();
    double getM();
    double getQ();
    int    getPDGID();
    TLorentzVector get4P();
};

// constructor and method definitions
inline leptons::leptons(double Pt, double Eta, double Phi, double M, double Q, int PDGID)
{
    pt=Pt;
    eta=Eta;
    phi=Phi;
    m=M;
    q=Q;
    pdgid=PDGID;
};

void leptons::setAll(double Pt, double Eta, double Phi, double M, double Q, int PDGID)
{
    pt=Pt;
    eta=Eta;
    phi=Phi;
    m=M;
    q=Q;
    pdgid=PDGID;
};

void leptons::setPt(double Pt){pt=Pt;}
void leptons::setEta(double Eta){eta=Eta;}
void leptons::setPhi(double Phi){phi=Phi;}
void leptons::setM(double M){m=M;}
void leptons::setQ(double Q){q=Q;}
void leptons::setPDGID(int PDGID){pdgid=PDGID;}

double leptons::getPt(){return pt;}
double leptons::getEta(){return eta;}
double leptons::getPhi(){return phi;}
double leptons::getM(){return m;}
double leptons::getQ(){return q;}
int    leptons::getPDGID(){return  pdgid;}

TLorentzVector leptons::get4P()
{
    TLorentzVector dummy4P;
    dummy4P.SetPtEtaPhiM(pt, eta, phi, m);
    return dummy4P;
}

//class definitons
class candidate
{
    double pt, eta, phi, m;

    public:
    candidate() = default;

    candidate(double, double, double, double);
    candidate(jets, jets);
    candidate(leptons, leptons);
    candidate(candidate, jets);

    void setAll(double, double, double, double);
    void setPt(double);
    void setEta(double);
    void setPhi(double);
    void setM(double);


    double getPt();
    double getEta();
    double getPhi();
    double getM();

    TLorentzVector get4P();
};

// constructor and method definitions
inline candidate::candidate(double Pt, double Eta, double Phi, double M)
{
    pt   = Pt;
    eta  = Eta;
    phi  = Phi;
    m    = M;
};

inline candidate::candidate(jets daughter1, jets daughter2)
{
    TLorentzVector candidate_4P;
    candidate_4P = daughter1.get4P()+daughter2.get4P();

    pt  = candidate_4P.Pt();
    eta = candidate_4P.Eta();
    phi = candidate_4P.Phi();
    m   = candidate_4P.M();
};

inline candidate::candidate(leptons daughter1, leptons daughter2)
{
    TLorentzVector candidate_4P;
    candidate_4P = daughter1.get4P()+daughter2.get4P();

    pt  = candidate_4P.Pt();
    eta = candidate_4P.Eta();
    phi = candidate_4P.Phi();
    m   = candidate_4P.M();
};

inline candidate::candidate(candidate daughter1, jets daughter2)
{
    TLorentzVector candidate_4P;
    candidate_4P = daughter1.get4P()+daughter2.get4P();

    pt  = candidate_4P.Pt();
    eta = candidate_4P.Eta();
    phi = candidate_4P.Phi();
    m   = candidate_4P.M();
};


void candidate::setAll(double Pt, double Eta, double Phi, double M)
{
    pt   = Pt;
    eta  = Eta;
    phi  = Phi;
    m    = M;
};

void candidate::setPt(double Pt){pt=Pt;}
void candidate::setEta(double Eta){eta=Eta;}
void candidate::setPhi(double Phi){phi=Phi;}
void candidate::setM(double M){m=M;}

double candidate::getPt(){return pt;}
double candidate::getEta(){return eta;}
double candidate::getPhi(){return phi;}
double candidate::getM(){return m;}

TLorentzVector candidate::get4P()
{
    TLorentzVector dummy4P;
    dummy4P.SetPtEtaPhiM(pt, eta, phi, m);
    return dummy4P;
}
