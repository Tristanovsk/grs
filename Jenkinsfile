/*******************************
* Constantes
********************************/
// Création d'une variable "UUID" pour spécifier un ID unique "server-id" du serveur de test configuré au niveau de jFrog-CLI
def SERVERID = UUID.randomUUID().toString()

// CREDENTIALS
def artifactoryCredentials = "OBS2CO_ARTIFACTORY_CREDENTIALS"


pipeline {

    // Agent desactive car ils sont renseignes aux niveaux des stages selon le besoin
    agent none
    options {
        // Ne sauvegarde sur Jenkins que les 20 derniers builds
        buildDiscarder(logRotator(numToKeepStr: '20'))
        // Stop le build si il depasse 1H
        timeout(time: 1, unit: 'HOURS')
        // Desactive le checkout par defaut car effectue dans un stage plus bas
        skipDefaultCheckout()
        // Interdit d'avoir 2 execution en meme temps
        disableConcurrentBuilds()
    }

    environment {
        ARTI_TOKEN = credentials('OBS2CO_ARTIFACTORY_CREDENTIALS')
        ARTI_URL = "https://${artifactory_host}/artifactory"
        SONAR_TOKEN=credentials('OBS2CO_SONAR_CREDENTIALS')
        // Specifie le chemin contenant des fichiers pour l'utilisation de sonar scanner. Jenkins se basera sur le workspace comme base de chemin
        SONAR_SCANNER_OPTS="-Duser.home=.sonar-scanner"
        // Specifie le chemin de configuration de jfrog cli. Jenkins se basera sur le workspace comme base de chemin
        JFROG_CLI_HOME_DIR=".jfrog"
        // Permet de minimiser les logs inutiles ex : % d'avancement de l'upload ou download
        CI="true"
        // Construit l'URL du build Artifactory, Il faut encoder URL pour que ce soit bien reconnu par Artifactory
        JFROG_CLI_BUILD_NAME_URL = java.net.URLEncoder.encode(JFROG_CLI_BUILD_NAME, "UTF-8")
        ARTIFACTORY_BUILD_URL = "https://${artifactory_host}/artifactory/webapp/#/builds/${JFROG_CLI_BUILD_NAME_URL}/${JFROG_CLI_BUILD_NUMBER}"
        // Recuperation d'un credential Jenkins
        PROXY_TOKEN=credentials('obs2co-proxy')
        // Utilisation d'un compte projet pour le proxy CNES
        HTTP_PROXY = "http://${PROXY_TOKEN}@proxy-tech-web.cnes.fr:8060"
        HTTPS_PROXY = "http://${PROXY_TOKEN}@proxy-tech-web.cnes.fr:8060"
        NO_PROXY='cnes.fr'
        // Build name pour jfrog cli. Permet de retrouver les informations de build dans Artifactory. Cette cle doit etre unique
        // Utilise le nom du job Jenkins pour le retrouver plus simplement et assure une meilleure unicite de la cle
        JFROG_CLI_BUILD_NAME = "${env.JOB_NAME}".replace('%2F', ':')
        // Build number pour jfrog cli. Utilise le numero de build de Jenkins pour faciliter le lien entre le build Jenkins et celui d'Artifactory
        JFROG_CLI_BUILD_NUMBER = "${env.BUILD_NUMBER}"
        // Variable definissant toutes les variables d'environnement qui ne seront pas transmis à Artifactory. Il faut exclure les credentials
        // Liste non exhaustive
        JFROG_CLI_ENV_EXCLUDE = "*token*;*cred*;*proxy*;*secret*;*key*;*password*"
        sonarqube_host = sonarqube.cnes.fr
        artifactory_host = "https://artifactory.cnes.fr/"
        artifactoryUrl = "https://artifactory.cnes.fr/artifactory"
        artifactoryRegistryUrl = "artifactory.cnes.fr/docker/"
        gitlab_host = "https://gitlab.cnes.fr/"
    }

    // Declenchement automatique tous les samedis matins pour des besoins usine logicielle
    triggers { cron('H 6 * * 6') }
	
    parameters{
        booleanParam(name: 'DEBUG', defaultValue: false, description: 'Mode debug, ne nettoie pas le workspace à la fin')
        booleanParam(name: 'LAUNCH_XRAY', defaultValue: true, description: 'Lance les scans Xray dans le job')
    }


    stages {
        /** Ce job lance une analyse qualite via le module hadolint qui n'est pas disponible sur les agents Jenkins docker car ces agents ne doivent servir que pour les commandes docker.
         *  L'analyse qualite est donc realise sur un noeud HPC. Les etapes sont donc parallelises pour gagner du temps
        **/
        stage('Parallel stage') {
            parallel {
                stage('Build docker') {
                    // Execution de la pipeline sur un agent docker
                    agent { label 'docker' }

                    stages {
                        stage ('clean workspace') {
                            steps {
                                // Actions faite au début au cas ou le précédent build avait l'option DEBUG d'activé
                                cleanWs()
                                // Nettoyage des éléments docker tel que les images, les containers, les réseaux et les volumes.
                                sh """
                                docker system prune -af
                                """
                            }
                        }

                        stage('checkout') {
                            steps {
                                // Clone Git renseigne lors de la configuration du job Jenkins au niveau de l'IHM
                                checkout scm
                                sh ls
                                sh echo "toto"
                            }
                        }

                        stage ("configuration Jfrog cli serveur") {
                            steps {
                                sh '''
                                # Configuration du serveur Artifactory
                                jfrog rt config '''+SERVERID+''' --url='''+ARTI_URL+''' --user=$ARTI_TOKEN_USR --access-token=$ARTI_TOKEN_PSW --enc-password --interactive=false

                                # Enregistrement des informations Git
                                jfrog rt build-add-git --server-id '''+SERVERID+'''

                                # Collecte des variables d'environnement pour ensuite les ajouter aux informations de construction:
                                jfrog rt bce
                                '''
                            }
                        }

                        stage('build') {
                            steps {
                                sh """
                                    export no_proxy=cnes.fr; curl -v -u '${username}:${token}' --insecure -O                                     
                                    docker login docker.pkg.github.com --username cecile.betmont --password ${token}
                                    docker pull docker.pkg.github.com/snap-contrib/docker-snap/snap:latest
                                    docker tag docker-snap/snap ${artifactoryRegistryUrl}/obs2co-docker-local/snap:latest
                                    jfrog rt docker-push --skip-login --server-id ${SERVERID} ${artifactoryRegistryUrl}/obs2co-docker-local/snap:latest
           
                                    #enregistrement des credentials proxy dans un fichier texte qui sera transmis à l'image docker
                                    echo http://${PROXY_TOKEN_USR}:${PROXY_TOKEN_PSW}@proxy-tech-web.cnes.fr:8060 > ./http_proxy.txt
                                    echo http://${PROXY_TOKEN_USR}:${PROXY_TOKEN_PSW}@proxy-tech-web.cnes.fr:8060 > ./https_proxy.txt

                                """

                                script {
                                    docker.withRegistry("https://${artifactory_host}/artifactory", 'obs2co-docker-local') {
                                        sh """
                                            mkdir -p certs
                                        #copie des certificats de l'agent docker dans le dossier certs/ pour ensuite les intégrer dans l'image Docker
                                            cp /etc/pki/ca-trust/source/anchors/AC*.crt certs/
                                        #transmission des credentials proxy à l'image en passant par le système de secrets
                                            DOCKER_BUILDKIT=1 docker build -t ${artifactory_host}/obs2co-docker-local/grs:latest --no-cache \
                                            --build-arg IMAGE_SOURCE=${artifactory_host}/obs2co-docker-local/snap \
                                            --build-arg no_proxy=cnes.fr \
                                            --secret id=proxy_http_cnes,src=http_proxy.txt \
                                            --secret id=proxy_https_cnes,src=https_proxy.txt \
                                            .
                                        """
                                    }
                                }
                            }
                        }

                        stage('livraison') {
                            steps {
                                script {
                                    docker.withRegistry("https://${artifactory_host}/artifactory", 'obs2co-docker-local') {
                                        sh  """
                                        # Publie sur Artifactory
                                        jfrog rt docker-push --skip-login --server-id ${SERVERID} ${artifactory_host}/obs2co-docker-local/example-docker-icode:latest testci-docker

                                        # Publication de l'objet build-info dans Artifactory. La variable BUILD_URL est une variable defini par Jenkins.
                                        jfrog rt bp --server-id ${SERVERID} --build-url ${BUILD_URL}
                                        """
                                    }

                                    // Permet d'afficher une icone contenant l'URL du build Artifactory directement dans l'historique des builds Jenkins
                                    currentBuild.description = "<a href='${ARTIFACTORY_BUILD_URL}'><img src='/plugin/artifactory/images/artifactory-icon.png' alt='[Artifactory]' title='Artifactory Build Info' width='16' height='16'></a>"
                                }
                            }
                        }                        
                        
                        stage('analyse Xray') {
                            steps {
                                script {
                                    if (params.LAUNCH_XRAY) {
                                        // Analyse Xray de l'artefact. Ne fait pas echouer le build si le scan echoue grace a la commande --fail
                                        sh "jfrog rt bs --server-id ${SERVERID} --fail=false"
                                    }                                    
                                }
                            }
                        } 
                    }
                    

                    post { 
                        cleanup {
                            script {
                                if (!params.DEBUG) {
                                    // Nettoyage des éléments docker tel que les images, les containers, les réseaux et les volumes.
                                    sh """
                                    docker system prune -af
                                    """
                                    // Nettoyage du workspace ou est spécifié le serveur d'authentification "server-id Jfrog-CLI" au niveau de l'agent jenkins 
                                    // pour éviter que d'autres projets l'utilisent "Bonne pratique de sécurité"
                                    cleanWs()
                                }
                            }
                        }
                    }
                }

                stage ('Qualité') {
                    // Execution de la pipeline sur un agent docker
                    agent { label 'hpc' }
                    stages {
                        stage ('clean workspace') {
                            steps {
                                // Actions faite au début au cas ou le précédent build avait l'option DEBUG d'activé
                                cleanWs()
                            }
                        }

                        stage('checkout') {
                            steps {
                                // Clone Git renseigne lors de la configuration du job Jenkins au niveau de l'IHM
                                checkout scm
                            }
                        }

                        stage('analyse Sonarqube') {
                            steps {
                                withSonarQubeEnv(sonarqube_host) {
                                    script {
                                        def targetBranch = ""

                                        // Si le nom de branche commence par feature/, hotfix/ ou release/, on analyse une short-lived par rapport à master. 
                                        // Ceci car notre workflow est un feature branch flow
                                        if (BRANCH_NAME ==~ /^(feature|hotfix|release)\/.*$/) {
                                            targetBranch = "-Dsonar.branch.target=${sonarqube_main_branch}"
                                        } else if ((BRANCH_NAME != sonarqube_main_branch)) {
                                            // Si rien ne correspond, c'est que le nom de branche ne respecte pas le pattern de votre projet
                                            // Ou que le nom de la branche n'est pas spécifié dans le job
                                            error 'Current branch (${BRANCH_NAME}) does not match any pattern or is empty.'
                                        }
                                                                
                                        // Le premier parametre est la branche courante, le second la long-live branche
                                        def properties = " -Dsonar.branch.name=${env.BRANCH_NAME} ${targetBranch}"

                                        sh '''
                                        # Charge le module hadolint
                                        module load hadolint

                                        # Analyse Hadolint avec un rapport au format Checkstyle pour intégration dans SonarQube.
                                        # Nous devons faire un double pipe true car en cas d'erreur remonte dans le rapport, la commande hadolint remonte \$? = 1 et nous ne voulons pas arreter la pipeline pour cela
                                        hadolint -f checkstyle Dockerfile > hadolint-report.xml || true
                                        
                                        # Import du rapport dans SonarQube
                                        sonar-scanner -Dsonar.login=$SONAR_TOKEN '''+properties+'''
                                        '''
                                    }
                                }
                            }
                        }

                        stage("quality Gate") {
                            steps {
                                timeout(time: 1, unit: 'HOURS') {
                                    // Rend le build instable si le quality Gate de l'analyse SonarQube echoue. Ceci dans le but de continuer le job pour demontrer son fonctionnement complet
                                    catchError(stageResult: "FAILURE", buildResult: "UNSTABLE") {
                                        script {
                                            withSonarQubeEnv(sonarqube_host) {
                                                def ceTask 
                                                def ceTaskUrl = sh(returnStdout: true, script: 'ce=`grep ceTaskUrl= .scannerwork/report-task.txt`; echo ${ce:10}').trim()
                                                timeout(time: 5, unit: 'MINUTES') {
                                                    waitUntil {
                                                        ceTask = sh(returnStdout: true, script: 'curl -u $SONAR_TOKEN ' + ceTaskUrl).trim()
                                                        echo ceTask.toString()
                                                        return ceTask.contains("\"status\":\"SUCCESS\"")
                                                    }
                                                }
                                                def analysisId = sh(returnStdout: true, script: "echo ${ceTask} | grep -m1 -oP 'analysisId:\\s*\\K[^}]+'").trim()
                                                def qualityGateUrl = "https://${sonarqube_host}" + "/api/qualitygates/project_status?analysisId=" + analysisId
                                                def qualityGate = sh(returnStdout: true, script: 'curl -u $SONAR_TOKEN '+ qualityGateUrl).trim()
                                                echo qualityGate.toString()
                                                if (qualityGate.contains("ERROR")) {
                                                    throw new Exception("Quality Gate failed")
                                                } else {
                                                    echo "Quality Gate success"
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    post { 
                        cleanup {
                            script {
                                if (!params.DEBUG) {
                                    // Nettoyage du workspace ou est spécifié le serveur d'authentification "server-id Jfrog-CLI" au niveau de l'agent jenkins 
                                    // pour éviter que d'autres projets l'utilisent "Bonne pratique de sécurité"
                                    cleanWs()
                                }
                            }
                        }
                    }
                }
            }
        }     
    }

    post { 
        failure {
            script {
                if (!params.DEBUG) {
                    // Modifier l'adresse mail en fonction de la votre
                    mail to: 'cecile.betmont@thalesgroup.com',
                        subject: "Failed Pipeline: ${currentBuild.fullDisplayName}",
                        body: "Something is wrong with ${env.BUILD_URL}"
                }
                echo "I failed :( "
            }
        }
    }
}
