using CSV

using DataFrames

data = CSV.read("data.csv", DataFrame)
print(data)

function longueur(a,b)
  return sqrt((data[a,1]-data[b,1])^2+(data[a,2]-data[b,2])^2)
end

 longueur(1,2)

function vecteur_unitaire(a,b)
    n = zeros(Float64, 2)
    l=longueur(a,b)
    n[1] = (1/l)*(data[b,1] - data[a,1])
    n[2] = (1/l)*(data[b,2] - data[a,2])
    return n
end

vecteur_unitaire(1,2)

function n4(l1,l2)
    N =  zeros(Float64, 4,1)
    v= vecteur_unitaire(l1,l2)
    s= vecteur_unitaire(l2,l1)
    N[1,1]= v[1]
    N[2,1]= v[2]
    N[3,1]= s[1]
    N[4,1]= s[2]
    return N
end
    

using LinearAlgebra

function reg(A,E,n1,n2)
    regi=zeros(Float64, 4,4)
    s=zeros(Float64, 4,1)
    s=n4(n1,n2)
    l =longueur(n1,n2)
    regi=((E*A)/l)*(s.*transpose(s))
    return regi
end

function DDL()
    compt=0
    s =1.0
    x=data[:,9]
    for p in x
        if p==0
            compt+=2
        end
        if p==1 || p==2
            compt+=1
        end
    end
    #les vecteurs ddl
    vecteur=zeros(Float64,length(x),2)
    for i in 1:(length(x))
        k=x[i]
        if k==0
            vecteur[i,1]=s
            s=s+1
            vecteur[i,2]=s
            s=s+1
         end
        if k==1
            vecteur[i,1]=0
            vecteur[i,2]=s
            s=s+1
        end
        if k==2
            vecteur[i,1]=s
            s=s+1
            vecteur[i,2]=0
        end
        if k==3
            vecteur[i,1]=0
            vecteur[i,2]=0
        end
    end
    return vecteur
end



function ddlg(a , b,s)
    dd_ = zeros(Float64,4)
    for i in 1:2
        if s[a,i]!=0 
            dd_[i]=s[a,i]
        end
        if s[b,i]!=0
            dd_[2+i]=s[b,i]
        end
    end
    return dd_
end

function fct(i,j,ch)
    indice =[-1,-1,-1,-1,-1,-1,-1,-1,]
    id=[ ]
    f = ddlg(i ,j,DDL())
    for z in 1:4
        if(f[z]==1)
            indice[1]=z
        end
        if(f[z]==2)
            indice[2]=z
        end
        if(f[z]==3)
            indice[3]=z
        end
        if(f[z]==4)
            indice[4]=z
        end
        if(f[z]==5)
            indice[5]=z
        end
        if(f[z]==6)
            indice[6]=z
        end
        if(f[z]==7)
            indice[7]=z
        end
        if(f[z]==8)
            indice[8]=z
        end
    end
    v=1
    if ch[1] in f && ch[2] in f 
        for z in ch
            if indice[z] != -1.0
                if (indice[z] ∉ id)
                    append!(id, indice[z])
                    v=v+1
                end
            end
        end
    end
    if ch[1]==ch[2] && length(id)>0
         append!(id,id[1])
    end
    return id
end

liaison=[[1,2],[1,4],[1,5],[1,6],[2,3],[2,4],[2,5],[2,6],[3,4],[3,5],[3,6],[4,5],[5,6]]

#matrice d'assemblage 
kl = zeros(Float64,8,8)
for k in 1:8 
    for b in k:8 
        for i in liaison 
            f=fct(i[1],i[2],[k,b])
            t=length(f)
            if t>1 
                m =reg(100*10^(-6),200000*10^6,i[1],i[2])
                kl[k,b]=kl[k,b]+m[f[1],f[2]] 
            end
        end
     end
end
for i in 1:8
    for j in (i+1):8
        kl[j,i]= kl[i,j]
    end
end
print("matrice d'assemblage") 
print("\n")
for i in 1:8
    print("[ ")
    for j in 1:8
        print(kl[i,j],"    ")
    end
    print("]")
    print("\n")
    print("\n")
end


using LinearAlgebra
b=zeros(Float64,8)
b[2]=-2000

print("les déplacements inconnus")
solution=kl\b

function energie(a,b,x)
    u = [ ]
    s=ddlg(a , b,DDL())
    d= Dict(0=>0,1=>solution[1],2=>solution[2],3=>solution[3],4=>solution[4],5=>solution[5],6=>solution[6],7=>solution[7],8=>solution[8])
    for i in s
        append!(u,d[i])
    end
    Edef = u'*reg(x,200000*10^6,a,b)*u
    return (1/2)*Edef
end

s=0
for i in liaison
    s=s+energie(i[1],i[2],10^(-4))
    print("energie de deformation de la barre liaison entre ",i[1],i[2],": ",energie(i[1],i[2],10^(-4)),"\n")
end
print("energie globale de la structure : ",s)

