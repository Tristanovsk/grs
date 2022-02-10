/*******************************
* Constantes
********************************/
// Création d'une variable "UUID" pour spécifier un ID unique "server-id" du serveur de test configuré au niveau de jFrog-CLI
def SERVERID = UUID.randomUUID().toString()

// CREDENTIALS
def artifactoryCredentials = "OBS2CO_ARTIFACTORY_TOKEN"


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
        ARTI_TOKEN = credentials('OBS2CO_ARTIFACTORY_TOKEN')
        DOCKER_TOKEN = credentials('DOCKER_TOKEN')
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
        sonarqube_host = "sonarqube.cnes.fr"
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
                                sh 'ls'
                                sh 'echo "toto"'
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
                                    export http_proxy='http://${PROXY_TOKEN}@proxy-tech-web.cnes.fr:8060'                             
                                    export https_proxy='http://${PROXY_TOKEN}@proxy-tech-web.cnes.fr:8060'                                   
                                    #docker pull artifactory.cnes.fr/obs2co-docker/snap-contrib/docker-snap/snap:latest
                                    
                                    #enregistrement des credentials proxy dans un fichier texte qui sera transmis à l'image docker
                                    echo http://${PROXY_TOKEN_USR}:${PROXY_TOKEN_PSW}@proxy-tech-web.cnes.fr:8060 > ./http_proxy.txt
                                    echo http://${PROXY_TOKEN_USR}:${PROXY_TOKEN_PSW}@proxy-tech-web.cnes.fr:8060 > ./https_proxy.txt

                                """
                                script {
                                    docker.withRegistry("${artifactory_host}/artifactory", 'OBS2CO_ARTIFACTORY_TOKEN') {
                                        sh """
                                            mkdir -p certs
                                        #copie des certificats de l'agent docker dans le dossier certs/ pour ensuite les intégrer dans l'image Docker
                                            cp /etc/pki/ca-trust/source/anchors/AC*.crt certs/
                                        #transmission des credentials proxy à l'image en passant par le système de secrets
                                            DOCKER_BUILDKIT=1 docker build -t artifactory.cnes.fr/obs2co-docker/grs:latest --no-cache \
                                            --build-arg IMAGE_SOURCE=artifactory.cnes.fr/obs2co-docker/snap-contrib/docker-snap \
                                            --build-arg no_proxy=cnes.fr \
                                            --secret id=proxy_http_cnes,src=http_proxy.txt \
                                            --secret id=proxy_https_cnes,src=https_proxy.txt \
                                            --build-arg HTTP_PROXY='http://${PROXY_TOKEN}@proxy-tech-web.cnes.fr:8060 \
                                            --build-arg HTTP_PROXY='http://${PROXY_TOKEN}@proxy-tech-web.cnes.fr:8060 \
                                            --build-arg USER=$ARTI_TOKEN_USR \
                                            --build-arg PSW=$ARTI_TOKEN_PSW \
                                            .
                                        """
                                    }
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
