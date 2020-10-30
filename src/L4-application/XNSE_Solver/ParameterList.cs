using BoSSS.Foundation;
using Microsoft.SqlServer.Server;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver
{
    interface IParameter
    {
        IList<string> Names { get; }

        DelParameterFactory Factory { get; }

        DelPartialParameterUpdate Update { get; }
    }

    class ParameterList
    {
        List<IParameter> parameters;

        public ParameterList()
        {
            parameters = new List<IParameter>(10);
        }

        public void AddParameter(IParameter parameter)
        {
            parameters.Add(parameter);
        }

        public void AddParameter(DelParameterFactory factory, DelPartialParameterUpdate update, IList<string> names)
        {
            Parameter parameter = new Parameter
            {
                Names = names,
                Factory = factory,
                Update = update
            };
            AddParameter(parameter);
        }

        class Parameter :IParameter
        {
            public IList<string> Names { get; set; }

            public DelParameterFactory Factory { get; set; }

            public DelPartialParameterUpdate Update { get; set; }
        }

        public ICollection<DelParameterFactory> Factories(IList<string> names)
        {
            LinkedList<string> nameList = new LinkedList<string>(names);
            LinkedList<DelParameterFactory> parameterFactories = new LinkedList<DelParameterFactory>();
            //Find parameters and remove all found parameters from list;

            while(nameList.Count > 0)
            {
                string name = nameList.First.Value;
                nameList.RemoveFirst();
                //Find currentName
                for(int i = 0; i < parameters.Count; ++i)
                {
                    IParameter parameter = parameters[i];
                    if (parameter.Names.Contains(name))
                    {
                        if(parameter.Factory != null)
                        {
                            parameterFactories.AddLast(parameter.Factory);
                        }
                        foreach(string otherParamName in parameter.Names)
                        {
                            nameList.Remove(otherParamName);
                        }
                        break;
                    }
                }
            }
            return parameterFactories;
            //AddTo Array

            
        }

        public ICollection<DelPartialParameterUpdate> ParameterUpdates(IList<string> names)
        {
            LinkedList<string> nameList = new LinkedList<string>(names);
            LinkedList<DelPartialParameterUpdate> parameterUpdates = new LinkedList<DelPartialParameterUpdate>();
            
            //Find parameters and remove all found parameters from list;
            while (nameList.Count > 0)
            {
                string name = nameList.First.Value;
                nameList.RemoveFirst();
                //Find currentName
                for (int i = 0; i < parameters.Count; ++i)
                {
                    IParameter parameter = parameters[i];
                    if (parameter.Names.Contains(name))
                    {
                        if(parameter.Update != null)
                        {
                            parameterUpdates.AddLast(parameter.Update);
                        }
                        foreach (string otherParamName in parameter.Names)
                        {
                            nameList.Remove(otherParamName);
                        }
                        break;
                    }
                }
            }
            return parameterUpdates;
        }
    }
}
