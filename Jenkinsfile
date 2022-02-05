/*******************************
* Constantes
********************************/
// Création d'une variable "UUID" pour spécifier un ID unique "server-id" du serveur de test configuré au niveau de jFrog-CLI
def SERVERID = UUID.randomUUID().toString()

def artifactoryUrl = "https://artifactory.cnes.fr/artifactory"
def artifactoryRegistryUrl = "artifactory.cnes.fr/docker/"

// CREDENTIALS
def artifactoryCredentials = "OBS2CO_ARTIFACTORY_CREDENTIALS"

// PROJECT VARIABLES
def projectName = "WaterQuality/grs2"
def deliveryPath = "obs2co-docker-local"


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
        ARTI_TOKEN = credentials('obs2co-docker-local')
        ARTI_URL = "https://${artifactory_host}/artifactory"
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
        JFROG_CLI_BUILD_NAME = "${env.JOB_NAME}".replace('%2F', ':')
        // Build number pour jfrog cli. Utilise le numero de build de Jenkins pour faciliter le lien entre le build Jenkins et celui d'Artifactory
        JFROG_CLI_BUILD_NUMBER = "${env.BUILD_NUMBER}"
        // Variable definissant toutes les variables d'environnement qui ne seront pas transmis à Artifactory. Il faut exclure les credentials
        // Liste non exhaustive
        JFROG_CLI_ENV_EXCLUDE = "*token*;*cred*;*proxy*;*secret*;*key*;*password*"
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
           stage('Init') {
                steps {
                        echo "========================== Init step =========================="

                }
           }
           stage('Build docker') {
                    // Execution de la pipeline sur un agent docker
                    agent { label 'docker' }

                    stages {
                        stage ('clean workspace') {
                            steps {
                                // Actions faite au début au cas ou le précédent build avait l'option DEBUG d'activé
                                cleanWs()
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

                        stage('pull') {
                            steps {
                                 withCredentials([usernamePassword(credentialsId: artifactoryCredentials, usernameVariable: 'username', passwordVariable: 'token')]) {
                                 sh """
                                    export no_proxy=cnes.fr; curl -v -u '${username}:${token}' --insecure -O                                     
                                    docker login docker.pkg.github.com --username cecile.betmont --password ${token}
                                    docker pull docker.pkg.github.com/snap-contrib/docker-snap/snap:latest
                                    jfrog rt docker-push --skip-login --server-id ${SERVERID} ${artifactoryRegistryUrl}/obs2co-docker-local/snap:latest snap
           
                                    #enregistrement des credentials proxy dans un fichier texte qui sera transmis à l'image docker
                                    echo http://${PROXY_TOKEN_USR}:${PROXY_TOKEN_PSW}@proxy-tech-web.cnes.fr:8060 > ./http_proxy.txt
                                    echo http://${PROXY_TOKEN_USR}:${PROXY_TOKEN_PSW}@proxy-tech-web.cnes.fr:8060 > ./https_proxy.txt
                                """
                                 }
                                script {
                                    docker.withRegistry("https://${artifactory_host}/artifactory", 'obs2co-docker-local') {
                                        sh """
                                            mkdir -p certs
                                        #copie des certificats de l'agent docker dans le dossier certs/ pour ensuite les intégrer dans l'image Docker
                                            cp /etc/pki/ca-trust/source/anchors/AC*.crt certs/
                                        #transmission des credentials proxy à l'image en passant par le système de secrets
                                            DOCKER_BUILDKIT=1 docker build -t ${artifactoryRegistryUrl}/obs2co-docker-local/grs:latest --no-cache \
                                            --build-arg IMAGE_SOURCE=${artifactoryRegistryUrl}/obs2co-docker-local/ \
                                            --build-arg no_proxy=cnes.fr \
                                            --secret id=proxy_http_cnes,src=http_proxy.txt \
                                            --secret id=proxy_https_cnes,src=https_proxy.txt \
                                            .
                                        """
                                    }
                                }
                            }
                        }


                        stage('build') {
                            steps {
                                sh """
                                    echo "Download le binaire pour l'image docker"
                                    git clone git@gitlab.cnes.fr/waterquality/grs2.git
                                    echo "Extraction du repo git de grs"
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
                                            DOCKER_BUILDKIT=1 docker build -t ${artifactoryRegistryUrl}/obs2co-docker-local/grs:latest --no-cache \
                                            --build-arg IMAGE_SOURCE=${artifactoryRegistryUrl}/obs2co-docker-local/snap \
                                            --build-arg no_proxy=cnes.fr \
                                            --secret id=proxy_http_cnes,src=http_proxy.txt \
                                            --secret id=proxy_https_cnes,src=https_proxy.txt \
                                            .
                                        """
                                    }
                                }
                            }
                        }
                    

                    
            }
        }
    }
}     
    



